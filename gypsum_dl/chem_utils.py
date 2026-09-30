"""The module includes definitions to manipulate the molecules."""

import copy
from typing import TYPE_CHECKING

from rdkit import Chem

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


def pick_lowest_enrgy_mols(
    mol_lst: list["MyMol"],
    num: int,
    thoroughness: int,
    energy_cache: dict[str, float | None] | None = None,
) -> list["MyMol"]:
    """Pick molecules with low energies. If necessary, the definition also
       makes a conformer without minimization (so not too computationally
       expensive).

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

    # First, generate 3D structures. How many? num * thoroughness. mols_3d is
    # a list of Gypsum-DL MyMol.MyMol objects.
    mols_3d = utils.random_sample(mol_lst, num * thoroughness, "")

    # Now get the energies
    data = []
    for i, mol in enumerate(mols_3d):
        energy = first_conf_energy(mol, energy_cache)

        if energy is not None:
            data.append((energy, i))

    data.sort()

    # Now keep only best top few.
    data = data[:num]

    # Keep just the mols there.
    return [mols_3d[d[1]] for d in data]


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

            # Pick the lowest-energy molecules. Note that this ranks candidates
            # using a throwaway conformer, so the molecules themselves are
            # returned in the state they arrived in. The container carries the
            # cache of those throwaway energies, so a variant that survives
            # this step is not embedded again when the next step ranks it.
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
