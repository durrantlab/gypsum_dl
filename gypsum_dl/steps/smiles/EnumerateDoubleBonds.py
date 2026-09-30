"""Module for enumerating unspecified double bonds (cis vs. trans)."""

import __future__

import copy
import itertools
import math
import random

import gypsum_dl.MolObjectHandling as MOH
import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl import MyMol, chem_utils, utils

try:
    from rdkit import Chem
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")


def enumerate_double_bonds(
    contnrs,
    max_variants_per_compound,
    thoroughness,
    num_procs,
    job_manager,
    parallelizer_obj,
):
    """Enumerates all possible cis-trans isomers. If the stereochemistry of a
       double bond is specified, it is not varied. All unspecified double bonds
       are varied.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: A list.
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
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param job_manager: The multithred mode to use.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    # No need to continue if none are requested.
    if max_variants_per_compound == 0:
        return

    utils.log("Enumerating all possible cis-trans isomers for all molecules...")

    # Group the molecule containers so they can be passed to the parallelizer.
    params = []
    for contnr in contnrs:
        params.extend(
            (mol, max_variants_per_compound, thoroughness) for mol in contnr.mols
        )
    params = tuple(params)

    # Ruin it through the parallelizer.
    tmp = []
    if parallelizer_obj is None:
        # MultiThreading with one processor is the same dispatch serial mode
        # uses, and it draws a seed per job; calling run_one directly left
        # this path drawing from whatever generator state happened to be in
        # place, so --random_seed did not reach it.
        tmp = Parallelizer.MultiThreading(params, 1, parallel_get_double_bonded)
    else:
        tmp = parallelizer_obj.run(
            params, parallel_get_double_bonded, num_procs, job_manager
        )

    # Remove Nones (failed molecules)
    clean = Parallelizer.strip_none(tmp)

    # Flatten the data into a single list.
    flat = Parallelizer.flatten_list(clean)

    # Containers that produced no isomers keep their existing structures.
    utils.carry_over_unrepresented(
        contnrs,
        flat,
        "double-bond variant",
        "(WARNING: Unable to generate double-bond variant)",
    )

    # No dedup over flat here: it mixes containers, so a global pass would
    # strip a second container's variants whenever two inputs share a SMILES.
    # bst_for_each_contnr_no_opt already dedups within each container.

    # Keep only the top few compound variants in each container, to prevent a
    # combinatorial explosion.
    chem_utils.bst_for_each_contnr_no_opt(
        contnrs,
        flat,
        max_variants_per_compound,
        thoroughness,
        variant_desc="cis-trans isomers",
    )


def sample_bond_dir_configs(num_bonds: int, cap: int) -> list[tuple[bool, ...]]:
    """Choose which up/down direction assignments to try for a set of bonds.

    Enumerating every assignment costs 2**num_bonds, and the number of bonds
    involved grows about four times faster than the number of double bonds
    being varied, so the full product can cost far more time and memory than
    the requested number of variants could ever justify. Sampling bounds that
    work while still reaching a good spread of the cis/trans forms the bonds
    can produce.

    Args:
        num_bonds: How many bonds will have their direction varied.
        cap: The largest number of assignments to return.

    Returns:
        A list of tuples of per-bond up (True) or down (False) flags: the full
            product when it fits within the cap, otherwise that many distinct
            assignments sampled from it.
    """

    space_size = 2**num_bonds
    if space_size <= cap:
        return list(itertools.product([True, False], repeat=num_bonds))

    if num_bonds < 63:
        # random.sample indexes into the population, so it needs the
        # population's length to fit in a C ssize_t.
        masks = random.sample(range(space_size), cap)
    else:
        # The space is so much larger than the cap here that repeated draws
        # essentially never collide, but track them so the sample stays
        # distinct regardless.
        masks = []
        seen = set()
        while len(masks) < cap:
            mask = random.getrandbits(num_bonds)
            if mask not in seen:
                seen.add(mask)
                masks.append(mask)

    return [tuple(bool(mask >> i & 1) for i in range(num_bonds)) for mask in masks]


def _could_carry_stereochemistry(mol_with_hs: "Chem.Mol", bond_idx: int) -> bool:
    """Report whether a double bond has distinguishable substituents on each end.

    The candidate list is built from `GetStereo() is STEREONONE`, which is also
    true of bonds that can never be cis/trans (carbonyls, terminal alkenes).
    Counting those toward the sampling budget let a genuinely stereogenic bond
    be discarded before any isomer was generated.

    Args:
        mol_with_hs: Molecule with explicit hydrogens, so terminal heavy atoms
            are distinguishable from substituted ones.
        bond_idx: Index of the double bond to test.

    Returns:
        False when either end has no other neighbor (C=O) or carries two
            hydrogens (=CH2), otherwise True.
    """
    bond = mol_with_hs.GetBondWithIdx(bond_idx)
    ends = (
        (bond.GetBeginAtom(), bond.GetEndAtomIdx()),
        (bond.GetEndAtom(), bond.GetBeginAtomIdx()),
    )
    for atom, partner_idx in ends:
        others = [n for n in atom.GetNeighbors() if n.GetIdx() != partner_idx]
        # A lone hydrogen still counts (an N-H imine is E/Z), so only an end
        # whose two substituents are both hydrogen is ruled out.
        if not others or sum(1 for n in others if n.GetAtomicNum() == 1) >= 2:
            return False
    return True


def parallel_get_double_bonded(mol, max_variants_per_compound, thoroughness):
    """A parallelizable function for enumerating double bonds.

    :param mol: The molecule with a potentially unspecified double bond.
    :type mol: MyMol.MyMol
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
    :return: [description]
    :rtype: [type]
    """

    # For this to work, you need to have explicit hydrogens in place. Use a
    # local copy so the molecule stored in the container is not mutated (which
    # would make results depend on the job manager and desync the cached SMILES).
    rdkit_mol_with_hs = Chem.AddHs(mol.rdkit_mol)

    # Get all double bonds that don't have defined stereochemistry. Note that
    # these are the bond indexes, not the atom indexes. AddHs preserves the
    # existing bond indexes, so these remain valid on rdkit_mol_with_hs.
    unasignd_dbl_bnd_idxs = mol.get_double_bonds_without_stereochemistry()

    if len(unasignd_dbl_bnd_idxs) == 0:
        # There are no unassigned double bonds, so move on.
        return [mol]

    # Throw out any bond that is in a small ring.
    unasignd_dbl_bnd_idxs = [
        i
        for i in unasignd_dbl_bnd_idxs
        if not rdkit_mol_with_hs.GetBondWithIdx(i).IsInRingSize(3)
    ]
    unasignd_dbl_bnd_idxs = [
        i
        for i in unasignd_dbl_bnd_idxs
        if not rdkit_mol_with_hs.GetBondWithIdx(i).IsInRingSize(4)
    ]
    unasignd_dbl_bnd_idxs = [
        i
        for i in unasignd_dbl_bnd_idxs
        if not rdkit_mol_with_hs.GetBondWithIdx(i).IsInRingSize(5)
    ]
    unasignd_dbl_bnd_idxs = [
        i
        for i in unasignd_dbl_bnd_idxs
        if not rdkit_mol_with_hs.GetBondWithIdx(i).IsInRingSize(6)
    ]
    unasignd_dbl_bnd_idxs = [
        i
        for i in unasignd_dbl_bnd_idxs
        if not rdkit_mol_with_hs.GetBondWithIdx(i).IsInRingSize(7)
    ]

    # Drop the bonds that cannot be cis/trans before the budget below is spent
    # on them. The loop further down already skips them, but by then they have
    # displaced stereogenic bonds from the random selection.
    unasignd_dbl_bnd_idxs = [
        i
        for i in unasignd_dbl_bnd_idxs
        if _could_carry_stereochemistry(rdkit_mol_with_hs, i)
    ]

    if len(unasignd_dbl_bnd_idxs) == 0:
        return [mol]

    # Previously, I fully enumerated all double bonds. When there are many
    # such bonds, that leads to a combinatorial explosion that causes problems
    # in terms of speed and memory. Now, enumerate only enough bonds to make
    # sure you generate at least thoroughness * max_variants_per_compound.
    unasignd_dbl_bnd_idxs_orig_count = len(unasignd_dbl_bnd_idxs)
    num_bonds_to_keep = max(
        1, int(math.ceil(math.log(thoroughness * max_variants_per_compound, 2)))
    )
    random.shuffle(unasignd_dbl_bnd_idxs)
    unasignd_dbl_bnd_idxs = sorted(unasignd_dbl_bnd_idxs[:num_bonds_to_keep])

    # Get a list of all the single bonds that come off each double-bond atom.
    all_sngl_bnd_idxs = set([])
    dbl_bnd_count = 0
    for dbl_bnd_idx in unasignd_dbl_bnd_idxs:
        bond = rdkit_mol_with_hs.GetBondWithIdx(dbl_bnd_idx)

        atom1 = bond.GetBeginAtom()
        atom1_bonds = atom1.GetBonds()
        if len(atom1_bonds) == 1:
            # The only bond is the one you already know about. So don't save.
            continue

        atom2 = bond.GetEndAtom()
        atom2_bonds = atom2.GetBonds()
        if len(atom2_bonds) == 1:
            # The only bond is the one you already know about. So don't save.
            continue

        dbl_bnd_count = dbl_bnd_count + 1

        # Suffice it to say, RDKit does not deal with cis-trans isomerization
        # in an intuitive way...
        idxs_of_other_bnds_frm_atm1 = [b.GetIdx() for b in atom1.GetBonds()]
        idxs_of_other_bnds_frm_atm1.remove(dbl_bnd_idx)

        idxs_of_other_bnds_frm_atm2 = [b.GetIdx() for b in atom2.GetBonds()]
        idxs_of_other_bnds_frm_atm2.remove(dbl_bnd_idx)

        all_sngl_bnd_idxs |= set(idxs_of_other_bnds_frm_atm1)
        all_sngl_bnd_idxs |= set(idxs_of_other_bnds_frm_atm2)

    # Now come up with up/down combinations for those bonds. Each retained
    # double bond contributes up to four single bonds, so the full product is
    # up to 2**(4 * dbl_bnd_count) even though num_bonds_to_keep is only
    # logarithmic in the variant budget. Because only a minority of the
    # direction assignments leave every double bond fully specified, sample
    # generously relative to the number of variants requested. The floor keeps
    # the cheap cases (one or two double bonds) fully enumerated.
    all_sngl_bnd_idxs = list(all_sngl_bnd_idxs)
    all_atom_config_options = sample_bond_dir_configs(
        len(all_sngl_bnd_idxs),
        max(1024, 64 * thoroughness * max_variants_per_compound),
    )

    # Let the user know.
    if dbl_bnd_count > 0:
        utils.log(
            "\t"
            + str(mol.smiles(True))
            + " has "
            # + str(dbl_bnd_count)
            + str(
                # Not exactly right, I think, because should be dbl_bnd_count, but ok.
                unasignd_dbl_bnd_idxs_orig_count
            )
            + " double bond(s) with unspecified stereochemistry."
        )

        # The bonds dropped above are never assigned a direction, so the 3D
        # embedder picks their geometry arbitrarily later on. Say so, rather
        # than letting the count above imply that all of them were enumerated.
        if unasignd_dbl_bnd_idxs_orig_count > len(unasignd_dbl_bnd_idxs):
            utils.log(
                "\t\tTo avoid a combinatorial explosion, only "
                + str(len(unasignd_dbl_bnd_idxs))
                + " of the "
                + str(unasignd_dbl_bnd_idxs_orig_count)
                + " double bonds with unspecified stereochemistry were varied. "
                + "The stereochemistry of the others was left unspecified."
            )

    # Go through and consider each of the retained combinations.
    smiles_to_consider = set([])
    # The bond directions set below only become cis/trans stereo through the
    # legacy AssignStereochemistry.
    with MOH.legacy_stereo_perception():
        for atom_config_options in all_atom_config_options:
            # Make a copy of the original RDKit molecule.
            a_rd_mol = copy.copy(rdkit_mol_with_hs)
            # a_rd_mol = Chem.MolFromSmiles(mol.smiles())

            for bond_idx, direc in zip(all_sngl_bnd_idxs, atom_config_options):
                # Always done with reference to the atom in the double bond.
                if direc:
                    a_rd_mol.GetBondWithIdx(bond_idx).SetBondDir(
                        Chem.BondDir.ENDUPRIGHT
                    )
                else:
                    a_rd_mol.GetBondWithIdx(bond_idx).SetBondDir(
                        Chem.BondDir.ENDDOWNRIGHT
                    )

            # Assign the StereoChemistry. Required to actually set it.
            a_rd_mol.ClearComputedProps()
            Chem.AssignStereochemistry(a_rd_mol, force=True)

            # Add to list of ones to consider. Canonicalize without the hydrogens
            # added above: MyMol parses with sanitize=False, so any [H] in the
            # SMILES stays in the graph and in can_smi, the key used to tell
            # variants apart. That makes a molecule spelled with explicit
            # hydrogens count as distinct from the same molecule spelled without
            # them, wasting a max_variants_per_compound slot. RemoveHs keeps the
            # hydrogens that define a double bond's stereochemistry, so a few
            # necessarily remain.
            try:
                smiles_to_consider.add(
                    Chem.MolToSmiles(
                        Chem.RemoveHs(a_rd_mol), isomericSmiles=True, canonical=True
                    )
                )
            except Exception:
                # Some molecules still give troubles. Unfortunate, but these are
                # rare cases. Let's just skip these. Example:
                # CN1C2=C(C=CC=C2)C(C)(C)[C]1=[C]=[CH]C3=CC(=C(O)C(=C3)I)I
                continue

    # Remove ones that don't have "/" or "\". These are not real enumerated
    # ones. Sort so downstream selection order does not depend on set iteration
    # order (which varies with PYTHONHASHSEED across processes).
    smiles_to_consider = sorted(s for s in smiles_to_consider if "/" in s or "\\" in s)

    # Get the maximum number of / + \ in any string.
    cnts = [s.count("/") + s.count("\\") for s in smiles_to_consider]

    if not cnts:
        # There are no appropriate double bonds. Move on...
        return [mol]

    max_cnts = max(cnts)

    # Only keep those with that same max count. The others have double bonds
    # that remain unspecified.
    smiles_to_consider = [
        s[0] for s in zip(smiles_to_consider, cnts) if s[1] == max_cnts
    ]
    results = []
    for smile_to_consider in smiles_to_consider:
        # Make a new MyMol.MyMol object with the specified smiles.
        new_mol = MyMol.MyMol(smile_to_consider)

        # Sometimes you get an error if there's a bad structure otherwise. Add
        # the new molecule to the list of results, if it does not have a bizarre
        # substructure. Ask through the accessor: for a MyMol built from a
        # SMILES string, can_smi is still unset at this point, so reading the
        # cached attribute directly can never detect a failed canonicalization.
        if (
            new_mol.smiles()
            not in (
                False,
                None,
            )
            and not new_mol.remove_bizarre_substruc()
        ):
            # MyMol.__init__ points orig_smi and orig_smi_deslt at the
            # variant's own SMILES, so without this the field that traces a
            # pose back to the library entry just repeats the variant.
            new_mol.inherit_contnr_props(mol)
            new_mol.genealogy = mol.genealogy[:]
            new_mol.genealogy.append(
                f"{new_mol.smiles(True)} (cis-trans isomerization)"
            )
            results.append(new_mol)

    # Return the results.
    return results
