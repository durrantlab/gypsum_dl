"""The module includes definitions to manipulate the molecules."""

import copy
import math
from collections import Counter
from collections.abc import Hashable
from typing import TYPE_CHECKING, TypedDict

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from gypsum_dl import utils

if TYPE_CHECKING:
    # Annotations only. MolContainer imports this module, so keeping the import
    # out of the runtime graph leaves chem_utils importable on its own.
    from gypsum_dl.MyMol import MyMol


def first_conf_energy(
    mol: "MyMol", energy_cache: dict[str, float | None] | None = None
) -> float | None:
    """Report a molecule's first-conformer energy without changing it.

    Ranking a candidate must not alter it. Generating the conformer on the
    molecule itself reprotonates its rdkit_mol and leaves 3D coordinates
    attached, so variants that survived a pruning step would end up in a
    different state than variants from containers that never needed pruning
    (visible, for instance, as real 3D coordinates in a 2d_output_only run).
    Measuring on a throwaway copy keeps the decision free of side effects.

    Embedding that throwaway is the expensive part, and it runs in the
    dispatching process rather than through the parallelizer, so it is the one
    part of the SMILES pipeline that never spreads across cores. The optional
    cache is what keeps it from being paid again for a structure that has
    already been measured: a variant that survives one pruning step is a
    candidate again at the next one. Caching also fixes a structure's ranking
    energy for the life of the cache, so the same structure cannot be ranked
    two different ways within one compound.

    Args:
        mol: The candidate molecule to measure.
        energy_cache: Canonical SMILES to previously measured energy, or None
            to measure without caching. Only the embedded probe energy is
            cached; an energy read off coordinates the molecule already carries
            describes that specific coordinate set, not the structure.

    Returns:
        The energy of the first conformer, or None if no conformer could be
            generated.
    """

    if mol.conformers:
        # Coordinates already exist (the 3D steps build them before pruning),
        # so there is nothing to generate and nothing to protect.
        return mol.conformers[0].energy

    # smiles() reports a failed canonicalization as None, which is no kind of
    # cache key: two unrelated molecules would share it.
    key = mol.smiles()
    cacheable = energy_cache is not None and isinstance(key, str)

    if cacheable and key in energy_cache:
        return energy_cache[key]

    probe = copy.deepcopy(mol)
    probe.make_first_3d_conf_no_min()

    energy = probe.conformers[0].energy if probe.conformers else None

    if cacheable:
        energy_cache[key] = energy

    return energy


def _comparable_group_key(mol: "MyMol") -> Hashable:
    """Identify which candidates have UFF energies that can be compared.

    UFF energies share a scale only between molecules with the same atoms, and
    RDKit's UFF has no electrostatic term, so protonation states (which differ
    in hydrogen count and charge) cannot be ranked against one another.
    Tautomers and stereoisomers keep the formula and net charge of the state
    they came from, so this key also groups each protonation state with its
    descendants. A zwitterion and its neutral form share a key; they have the
    same atoms, and the tautomer step can interconvert them anyway.

    Args:
        mol: A candidate variant.

    Returns:
        The molecular formula and net formal charge, or the molecule itself
            when its SMILES cannot be determined or parsed, so such a molecule
            is a group of its own.
    """
    smi = mol.smiles()
    parsed = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if parsed is None:
        return mol
    return (rdMolDescriptors.CalcMolFormula(parsed), Chem.GetFormalCharge(parsed))


def pick_lowest_enrgy_mols(
    mol_lst: list["MyMol"],
    num: int,
    thoroughness: int,
    energy_cache: dict[str, float | None] | None = None,
) -> list["MyMol"]:
    """Pick molecules with low energies. If necessary, the definition also
       makes a conformer without minimization (so not too computationally
       expensive).

       Energies are compared only within groups of candidates that have the
       same formula and net charge (see _comparable_group_key). Slots are
       filled round-robin across groups: each group's lowest-energy candidate
       first, then each group's second lowest, and so on. Groups are visited
       in order of their lowest energy, which only matters when there are more
       groups than slots.

    :param mol_lst: The list of MyMol.MyMol objects.
    :type mol_lst: list
    :param num: The number of the lowest-energy ones to keep.
    :type num: int
    :param thoroughness: How many molecules to generate per variant (molecule)
       retained, for evaluation. For example, perhaps you want to advance five
       molecules (max_variants_per_compound = 5). You could just generate five
       and advance them all. Or you could generate ten and advance the best
       five (so thoroughness = 2). Using thoroughness > 1 increases the
       computational expense, but it also increases the chances of finding good
       molecules.
    :type thoroughness: int
    :param energy_cache: Ranking energies already measured for this compound,
       keyed on canonical SMILES, or None to measure without caching. See
       first_conf_energy. Defaults to None.
    :type energy_cache: dict, optional
    :return: Returns a list of MyMol.MyMol, the best ones.
    :rtype: list
    """

    # Remove identical entries. First-seen order, because set() order tracks
    # PYTHONHASHSEED and this list is the input to the seeded sampling below.
    mol_lst = list(dict.fromkeys(mol_lst))

    # If the length of the mol_lst is less than num, just return them all.
    if len(mol_lst) <= num:
        return mol_lst

    # Group in first-seen order, so a seeded run samples the groups in a fixed
    # order.
    groups: dict[Hashable, list["MyMol"]] = {}
    for mol in mol_lst:
        groups.setdefault(_comparable_group_key(mol), []).append(mol)

    # Split the embedding budget (num * thoroughness) across the groups. A
    # single sample over the pooled list could leave a protonation state
    # unsampled, and so eliminate it before it was ever scored.
    per_group = thoroughness * math.ceil(num / len(groups))

    # Score a sample of each group, lowest energy first. Ties keep sample
    # order, and candidates that could not be embedded are dropped.
    ranked_groups: list[list["MyMol"]] = []
    group_best: list[float] = []
    for members in groups.values():
        sample = utils.random_sample(members, per_group, "")
        scored: list[tuple[float, int]] = []
        for i, mol in enumerate(sample):
            energy = first_conf_energy(mol, energy_cache)
            if energy is not None:
                scored.append((energy, i))
        if scored:
            scored.sort()
            ranked_groups.append([sample[i] for _, i in scored])
            group_best.append(scored[0][0])

    # Stable sort, so groups with equal best energies keep first-seen order.
    order = sorted(range(len(ranked_groups)), key=lambda g: group_best[g])

    return _fill_round_robin(ranked_groups, order, num)


def _fill_round_robin(
    ranked_groups: list[list["MyMol"]], order: list[int], num: int
) -> list["MyMol"]:
    """Take each group's best candidate in turn, then each group's next best.

    Shared by both selectors, so that whichever ranks the candidates, every
    protonation state gets a slot before any state gets a second one.

    Args:
        ranked_groups: Each group's candidates, best first.
        order: The order in which to visit the groups.
        num: The number of candidates to keep.

    Returns:
        Up to num candidates.
    """
    kept: list["MyMol"] = []
    depth = 0
    while len(kept) < num and any(depth < len(g) for g in ranked_groups):
        for g in order:
            if depth < len(ranked_groups[g]) and len(kept) < num:
                kept.append(ranked_groups[g][depth])
        depth += 1

    return kept


def constitution_smiles(smiles: object) -> str | None:
    """Canonical SMILES without stereo, so tautomers compare by bonding only.

    Public because the tautomer worker keys its imine-shift bookkeeping on it,
    and the two sides must normalize SMILES the same way to match.

    Args:
        smiles: A SMILES string, or the None that smiles() reports on failure.

    Returns:
        The non-isomeric canonical SMILES, or None if there is none.
    """
    parsed = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if parsed is None:
        return None
    return Chem.MolToSmiles(parsed, isomericSmiles=False)


def tautomer_score(mol: "MyMol") -> int | None:
    """Score a tautomer with RDKit's rules for choosing a canonical tautomer.

    The rules reward features of the favored form, chiefly aromatic rings and
    carbonyls, and penalize forms such as hydrogens on sulfur or phosphorus.
    Unlike a force-field energy, the score can be compared across tautomers,
    which differ in bonding; and it needs no embedding.

    Args:
        mol: A candidate tautomer.

    Returns:
        The score (higher is better), or None if the molecule's SMILES cannot
            be determined or parsed.
    """
    smi = mol.smiles()
    parsed = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if parsed is None:
        return None
    return rdMolStandardize.TautomerEnumerator.ScoreTautomer(parsed)


# RDKit's tautomer scorer awards an oxime bonus (C=N-OH) that outweighs the
# carbonyl bonus, so on its own it prefers a hydroximic acid, C(OH)=N-OH, to
# the hydroxamic acid, C(=O)N-OH, that the group is in practice (compound 667
# of RDKit's bundled NCI sample; vorinostat and most HDAC inhibitors carry the
# group). Either oxygen may be deprotonated, since ionization runs first.
_HYDROXIMIC_ACID = Chem.MolFromSmarts(
    "[$([OX2H1]),$([OX1-])][CX3]=[NX2][$([OX2H1]),$([OX1-])]"
)

# A nitronic (aci-nitro) acid, C=N+(O-)OH, is a minor tautomer of its nitro
# compound by many orders of magnitude, yet the scorer can rank it above the
# nitro form. Without this, nizatidine's nitroethenediamine (CHEMBL3183075)
# became a nitronic acid once its imine form was demoted. The nitronate anion
# has no OH and is not matched.
_NITRONIC_ACID = Chem.MolFromSmarts("[CX3]=[N+]([OX1-])[OX2H1]")

# A beta-imino carbonyl, C(=N)-CH-C(=O), is the minor tautomer of a
# beta-enaminone or enamino ester, whose C=C-C=O is a vinylogous amide and
# whose N-H is hydrogen-bonded to the carbonyl. The imine/enamine demotion
# (MyMol.needs_imine_shift) misses forms MolVS reaches through keto/enol shifts
# instead, as for the Gould-Jacobs intermediate diethyl
# (3-chloroanilinomethylene)malonate (compound 3748 of RDKit's bundled NCI
# sample); naming the form catches it whatever the route. Excluded, because
# the imine form is the real one there: hydrazones and oximes (N-N or N-O on
# the imine nitrogen), and imine carbons that carry N, O, or S (amidines,
# imidates, and thioimidates such as 2-alkylthio-1-pyrrolines).
_BETA_IMINO_CARBONYL = Chem.MolFromSmarts(
    "[CX3;!$(C-[N,O,S])](=[NX2;!$(N-[N,O])])-[CX4;!H0]-[CX3]=[OX1]"
)

# The forms each overrated form should yield to. Either oxygen of a hydroxamic
# acid, and its nitrogen, may carry the charge after ionization.
_HYDROXAMIC_ACID = Chem.MolFromSmarts(
    "[CX3](=[OX1])[$([NX3]),$([NX2-])][$([OX2H1]),$([OX1-])]"
)
_NITRO = Chem.MolFromSmarts("[#6][N+](=[OX1])[OX1-]")
_BETA_ENAMINONE = Chem.MolFromSmarts("[NX3;!H0]-[CX3]=[CX3]-[CX3]=[OX1]")

# Each overrated form ranks below a candidate that has the form it should
# yield to, and only then. Demoting it below everything instead handed the slot
# to whatever else MolVS had made: compounds 1233 and 1367 of the NCI sample
# (2-(arylimino)methyl cyclohexanones) became ring enols, and compound 1905
# (2-thiouracil drawn as its 5,6-dihydro imine) became an iminol, because no
# enaminone was among their candidates to take the place of the demoted form.
_OVERRATED_FORMS: tuple[tuple[Chem.Mol, Chem.Mol], ...] = (
    (_HYDROXIMIC_ACID, _HYDROXAMIC_ACID),
    (_NITRONIC_ACID, _NITRO),
    (_BETA_IMINO_CARBONYL, _BETA_ENAMINONE),
)


class RuleHits(TypedDict):
    """Which of the _OVERRATED_FORMS pairs a candidate matches, by index."""

    overrated: frozenset[int]
    preferred: frozenset[int]


def overrated_form_hits(mol: "MyMol") -> RuleHits:
    """Find which overrated forms, and which preferred partners, a tautomer has.

    Demotion depends on the other candidates (see _OVERRATED_FORMS), so this
    reports both sides and leaves the decision to the selector.

    Args:
        mol: A candidate tautomer.

    Returns:
        The indices of the pairs whose overrated form the tautomer contains,
            and of those whose preferred form it contains. Both are empty if
            its SMILES cannot be parsed.
    """
    smi = mol.smiles()
    parsed = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
    if parsed is None:
        return {"overrated": frozenset(), "preferred": frozenset()}
    return {
        "overrated": frozenset(
            i
            for i, (bad, _) in enumerate(_OVERRATED_FORMS)
            if parsed.HasSubstructMatch(bad)
        ),
        "preferred": frozenset(
            i
            for i, (_, good) in enumerate(_OVERRATED_FORMS)
            if parsed.HasSubstructMatch(good)
        ),
    }


def pick_highest_scoring_tautomers(
    mol_lst: list["MyMol"], num: int, parent_smiles: set[str | None]
) -> list["MyMol"]:
    """Keep the tautomers RDKit's scoring rules rank highest.

    The tautomer step used pick_lowest_enrgy_mols through 2.0.0. Its energy is
    UFF on an unminimized embedding, and force-field energies cannot rank
    tautomers at all, because tautomers differ in bonding. On the first 1,000
    compounds of RDKit's bundled NCI sample, its pick matched the curated input
    form about 60% of the time, and the pick changed with the seed;
    tautomer_score matched about 92%, and it is deterministic. Ties go first
    to the form the step started from (the curated input, protonated), then to
    canonical SMILES.

    Two corrections come before the score, because the score has no term for
    what decides these cases. First, a tautomer with a hydroximic acid,
    nitronic acid, or beta-imino carbonyl group ranks below the other
    candidates when one of them has the form it should yield to (a hydroxamic
    acid, a nitro group, or a beta-enaminone); see _OVERRATED_FORMS.
    Second, a tautomer MolVS reached only through an imine/enamine shift ranks
    below the rest (MyMol.needs_imine_shift, set in parallel_make_taut). Which
    side of that pair dominates depends on conjugation: a simple imine beats
    its NH enamine, but an enamine conjugated with a carbonyl, nitrile, nitro
    group, or aryl ring (enaminones, enamino esters, arylidene hydantoins,
    nizatidine) beats its imine. The score always prefers the imine, while
    curated inputs draw whichever form is real, so the input's form is kept.
    The demoted form still fills a later slot when more variants are allowed.
    Not handled: a C=C that migrates along a chain through keto/enol shifts
    (diethyl itaconate to citraconate), which the score cannot tell apart.

    Candidates are still grouped by protonation state and
    filled round-robin, as in pick_lowest_enrgy_mols, with the groups visited
    in first-seen order: scores, like energies, do not compare across
    protonation states.

    Args:
        mol_lst: The candidate tautomers for one container.
        num: The number to keep.
        parent_smiles: SMILES of the variants the tautomers were generated
            from, for the tie-break. None entries (failed canonicalizations)
            are ignored.

    Returns:
        Up to num tautomers, best first within each protonation state.
            Candidates that cannot be scored are dropped, unless every
            candidate fits within num, in which case they come last.
    """

    # First-seen order, as in pick_lowest_enrgy_mols.
    mol_lst = list(dict.fromkeys(mol_lst))

    parents = {constitution_smiles(smi) for smi in parent_smiles} - {None}

    # Ranked even when every candidate fits, since scoring costs no embedding:
    # the order is what later steps and the output list, and a demoted form
    # belongs after the forms that outrank it whether or not anything is cut.
    unscorable: list["MyMol"] = []
    scored: dict[Hashable, list[tuple[int, RuleHits, "MyMol"]]] = {}
    for mol in mol_lst:
        score = tautomer_score(mol)
        if score is None:
            unscorable.append(mol)
            continue
        scored.setdefault(_comparable_group_key(mol), []).append(
            (score, overrated_form_hits(mol), mol)
        )

    ranked_groups: list[list["MyMol"]] = []
    for members in scored.values():
        # Partners count within the protonation group only, the same scope
        # the ranking itself has, and only on another candidate: a form cannot
        # yield to itself.
        partner_counts = Counter(i for _, h, _ in members for i in h["preferred"])
        keyed = []
        for score, hits, mol in members:
            smi = mol.smiles()
            elsewhere = {
                i
                for i in hits["overrated"]
                if partner_counts[i] - (1 if i in hits["preferred"] else 0) > 0
            }
            keyed.append(
                (
                    1 if elsewhere else 0,
                    1 if mol.needs_imine_shift else 0,
                    -score,
                    0 if constitution_smiles(smi) in parents else 1,
                    smi,
                    mol,
                )
            )
        ranked_groups.append([e[5] for e in sorted(keyed, key=lambda e: e[:5])])
    kept = _fill_round_robin(ranked_groups, list(range(len(ranked_groups))), num)

    # When nothing has to be cut, nothing is: an unscorable candidate is kept
    # after the ranked ones, as pick_lowest_enrgy_mols keeps every candidate
    # when they all fit.
    if len(mol_lst) <= num:
        kept += unscorable
    return kept


def remove_highly_charged_molecules(mol_lst):
    """Remove molecules that are highly charged.

    :param mol_lst: The list of molecules to consider.
    :type mol_lst: list
    :return: A list of molecules that are not too charged.
    :rtype: list
    """

    # First, find the molecule that is closest to being neutral. Ties on
    # |charge| (the common -1/+1 pair) are broken by preferring the more
    # negative form: a positional tie-break made the signed reference charge,
    # and so the kept set, depend on the order the variants happened to arrive
    # in.
    charges = [Chem.GetFormalCharge(mol.rdkit_mol) for mol in mol_lst]
    charge_closest_to_neutral = min(charges, key=lambda c: (abs(c), c))

    # Now create a new mol list, where the charges deviation from the most
    # neutral by no more than 4. Note that this used to be 2, but I increased
    # it to 4 to accommodate ATP.
    new_mol_lst = []
    for i, charge in enumerate(charges):
        if abs(charge - charge_closest_to_neutral) <= 4:
            new_mol_lst.append(mol_lst[i])
        else:
            # str() because smiles() reports failure as None, and this runs in
            # the main process, where a TypeError would end the run.
            utils.log(
                "\tWARNING: Discarding highly charged form: "
                + str(mol_lst[i].smiles())
                + "."
            )

    return new_mol_lst


def bst_for_each_contnr_no_opt(
    contnrs,
    mol_lst,
    max_variants_per_compound,
    thoroughness,
    crry_ovr_frm_lst_step_if_no_fnd=True,
    variant_desc: str = "variants",
    rank_by_tautomer_score: bool = False,
) -> None:
    """Keep only the top few compound variants in each container, to prevent a
       combinatorial explosion. This is run periodically on the growing
       containers to keep them in check.

       mol_lst mixes every container's candidates, and everything here treats
       it one container at a time, keyed on contnr_idx. A caller that
       deduplicates, samples, or otherwise operates on the flat list as a
       whole breaks that partition (two inputs that share a SMILES then starve
       each other), so callers must leave cross-container work to this
       function.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param mol_lst: The list of MyMol.MyMol objects.
    :type mol_lst: list
    :param max_variants_per_compound: To control the combinatorial explosion,
       only this number of variants (molecules) will be advanced to the next
       step.
    :type max_variants_per_compound: int
    :param thoroughness: How many molecules to generate per variant (molecule)
       retained, for evaluation. For example, perhaps you want to advance five
       molecules (max_variants_per_compound = 5). You could just generate five
       and advance them all. Or you could generate ten and advance the best
       five (so thoroughness = 2). Using thoroughness > 1 increases the
       computational expense, but it also increases the chances of finding good
       molecules.
    :type thoroughness: int
    :param crry_ovr_frm_lst_step_if_no_fnd: If it can't find any low-energy
       conformers, determines whether to just keep the old ones. Defaults to
       True.
    :param crry_ovr_frm_lst_step_if_no_fnd: bool, optional
    :param variant_desc: What the calling step was producing (e.g.,
       "tautomers"), for the warning logged when a container ends up with none.
       Defaults to "variants".
    :type variant_desc: str, optional
    :param rank_by_tautomer_score: Rank by pick_highest_scoring_tautomers
       instead of pick_lowest_enrgy_mols; thoroughness is then unused. Only
       the tautomer step sets this. Defaults to False.
    :type rank_by_tautomer_score: bool, optional
    :raises Exception: If a candidate's contnr_idx matches none of contnrs, or
       two containers share a contnr_idx.
    """

    # A candidate whose contnr_idx names no container would otherwise be
    # dropped without a word, leaving its real container to fall back on the
    # previous step's variants. contnrs_by_idx also rejects duplicate indices,
    # which would make the grouping below ambiguous.
    contnr_by_idx = utils.contnrs_by_idx(contnrs)
    for mol in mol_lst:
        if mol is not None and mol.contnr_idx not in contnr_by_idx:
            utils.exception(
                f"A candidate variant ({mol.smiles()}) carries contnr_idx "
                f"{mol.contnr_idx!r}, which matches none of the molecule "
                "containers."
            )

    # Remove duplicate ligands from each container.
    for mol_cont in contnrs:
        mol_cont.remove_identical_mols_from_contnr()

    # Group the smiles by contnr_idx.
    data = utils.group_mols_by_container_index(mol_lst)

    # Go through each container.
    for contnr_idx, contnr in enumerate(contnrs):
        contnr_idx = contnr.contnr_idx
        none_generated = False

        # Pick just the lowest-energy conformers from the new candidates.
        # Possible a compound was eliminated early on, so doesn't exist.
        if contnr_idx in list(data.keys()):
            mols = data[contnr_idx]

            # Remove molecules with unusually high charges.
            mols = remove_highly_charged_molecules(mols)

            if rank_by_tautomer_score:
                # contnr.mols still holds the variants the tautomers came
                # from; the container is repopulated only below.
                mols = pick_highest_scoring_tautomers(
                    mols,
                    max_variants_per_compound,
                    {m.smiles() for m in contnr.mols},
                )
            else:
                # Pick the lowest-energy molecules. Note that this ranks
                # candidates using a throwaway conformer, so the molecules
                # themselves are returned in the state they arrived in. The
                # container carries the cache of those throwaway energies, so
                # a variant that survives this step is not embedded again when
                # the next step ranks it.
                mols = pick_lowest_enrgy_mols(
                    mols,
                    max_variants_per_compound,
                    thoroughness,
                    contnr.probe_energies,
                )

            if len(mols) > 0:
                # Now remove all previously determined mols for this
                # container.
                contnr.mols = []

                # Add in the lowest-energy conformers back to the container.
                for mol in mols:
                    contnr.add_mol(mol)
            else:
                none_generated = True
        else:
            none_generated = True

        # No candidates survived for this container. Most callers are SMILES
        # steps, where no conformers are involved, so the message names what
        # the step was producing instead.
        if none_generated:
            no_survivors = (
                "\tWARNING: No "
                + variant_desc
                + " remained for "
                + contnr.orig_smi_deslt
                + " ("
                + contnr.name
                + ")."
            )
            if crry_ovr_frm_lst_step_if_no_fnd:
                # Just use previous ones.
                utils.log(
                    no_survivors + " Keeping the variants from the previous step."
                )
            else:
                # Discard the compound.
                utils.log(no_survivors + " Discarding the compound.")
                contnr.mols = []


def uniq_mols_in_list(mol_lst):
    # You need to make new molecules to get it to work.
    # new_smiles = [m.smiles() for m in self.mols]
    # new_mols = [Chem.MolFromSmiles(smi) for smi in new_smiles]
    # new_can_smiles = [Chem.MolToSmiles(new_mol, isomericSmiles=True, canonical=True) for new_mol in new_mols]

    can_smiles_already_set = set([])
    uniq_mols = []
    for m in mol_lst:
        smi = m.smiles()
        if not isinstance(smi, str):
            # smiles() reports failure as None. Two molecules that both failed
            # to canonicalize have not been shown to be the same molecule, so
            # keeping only the first would discard distinct structures; let the
            # steps that check for a usable molecule decide their fate.
            uniq_mols.append(m)
            continue
        if smi not in can_smiles_already_set:
            uniq_mols.append(m)
        can_smiles_already_set.add(smi)

        # if not new_can_smile in can_smiles_already_set:
        #     # Never seen before
        #     can_smiles_already_set.add(new_can_smile)
        # else:
        #     # Seen before. Delete!
        #     self.mols[i] = None

    # while None in self.mols:
    #     self.mols.remove(None)

    return uniq_mols
