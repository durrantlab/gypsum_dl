"""Module for enumerating unspecified double bonds (cis vs. trans)."""

import __future__

import itertools
import copy
import random

import gypsum.Parallelizer as Parallelizer
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

def enumerate_double_bonds(contnrs, max_variants_per_compound, thoroughness, num_procs, multithread_mode, parallelizer_obj):
    """Enumerates all possible cis-trans isomers. If the stereochemistry of a
       double bond is specified, it is not varied. All unspecified double bonds
       are varied.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: A list.
    :param max_variants_per_compound: [description] JDD: Figure out.
    :type max_variants_per_compound: int
    :param thoroughness: [description] JDD: Figure out.
    :type thoroughness: int
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param multithread_mode: The multithred mode to use.
    :type multithread_mode: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    # No need to continue if none are requested.
    if max_variants_per_compound == 0:
        return

    Utils.log(
        "Enumerating all possible cis-trans isomers for all molecules..."
    )

    # Group the molecule containers so they can be passed to the parallelizer.
    params = []
    for contnr in contnrs:
        for mol in contnr.mols:
            params.append(tuple([mol, max_variants_per_compound]))
    params = tuple(params)

    # Ruin it through the parallelizer.
    tmp = parallelizer_obj.run(params, parallel_get_double_bonded, num_procs, multithread_mode)

    # Remove Nones (failed molecules)
    clean = Parallelizer.strip_none(tmp)

    # Flatten the data into a single list.
    flat = Parallelizer.flatten_list(clean)

    # Get the indexes of the ones that failed to generate.
    contnr_idxs_of_failed = Utils.fnd_contnrs_not_represntd(contnrs, flat)

    # Go through the missing ones and throw a message.
    for miss_indx in contnr_idxs_of_failed:
        Utils.log(
            "\tCould not generate valid double-bond variant for " +
            contnrs[miss_indx].orig_smi + " (" +
            contnrs[miss_indx].name + "), so using existing " +
            "(unprocessed) structures.")
        for mol in contnrs[miss_indx].mols:
            mol.genealogy.append("(WARNING: Unable to generate double-bond variant)")
            clean.append(mol)

    # Keep only the top few compound variants in each container, to prevent a
    # combinatorial explosion.
    ChemUtils.bst_for_each_contnr_no_opt(contnrs, flat, max_variants_per_compound, thoroughness)

def parallel_get_double_bonded(mol, max_variants_per_compound):
    """[summary]

    :param mol: The molecule with a potentially unspecified double bond.
    :type mol: MyMol.MyMol
    :param max_variants_per_compound: [description] JDD: Figure out.
    :type max_variants_per_compound: int
    :return: [description]
    :rtype: [type]
    """

    # Get all double bonds that don't have defined stereochemistry. Note that
    # these are the bond indexes, not the atom indexes.
    unasignd_dbl_bnd_idxs = mol.get_double_bonds_without_stereochemistry()

    # Throw out any bond that is in a small ring.
    unasignd_dbl_bnd_idxs = [i for i in unasignd_dbl_bnd_idxs
                if not mol.rdkit_mol.GetBondWithIdx(i).IsInRingSize(3)]
    unasignd_dbl_bnd_idxs = [i for i in unasignd_dbl_bnd_idxs
                if not mol.rdkit_mol.GetBondWithIdx(i).IsInRingSize(4)]
    unasignd_dbl_bnd_idxs = [i for i in unasignd_dbl_bnd_idxs
                if not mol.rdkit_mol.GetBondWithIdx(i).IsInRingSize(5)]
    unasignd_dbl_bnd_idxs = [i for i in unasignd_dbl_bnd_idxs
                if not mol.rdkit_mol.GetBondWithIdx(i).IsInRingSize(6)]
    unasignd_dbl_bnd_idxs = [i for i in unasignd_dbl_bnd_idxs
                if not mol.rdkit_mol.GetBondWithIdx(i).IsInRingSize(7)]

    # Throw out any bond that has an atom that only participate in that one
    # bond (terminal alkene).
    idx_of_othr_bnds_to_use = []
    bonds_already_used = set([])
    print(Chem.MolToSmiles(mol.rdkit_mol))
    for dbl_bnd_idx in unasignd_dbl_bnd_idxs:
        bond = mol.rdkit_mol.GetBondWithIdx(dbl_bnd_idx)

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

        # Suffice it to say, RDKit does not deal with cis-trans isomerization
        # in an intuitive way...
        idxs_of_other_bnds_frm_atm1 = [b.GetIdx() for b in atom1.GetBonds()]
        idxs_of_other_bnds_frm_atm1.remove(dbl_bnd_idx)
        idxs_of_other_bnds_frm_atm1 = list(
            set(idxs_of_other_bnds_frm_atm1) - bonds_already_used
        )
        idxs_of_one_bnds_frm_atm1 = idxs_of_other_bnds_frm_atm1[0]

        idxs_of_other_bnds_frm_atm2 = [b.GetIdx() for b in atom2.GetBonds()]
        idxs_of_other_bnds_frm_atm2.remove(dbl_bnd_idx)
        idxs_of_other_bnds_frm_atm2 = list(
            set(idxs_of_other_bnds_frm_atm2) - bonds_already_used
            - set([idxs_of_one_bnds_frm_atm1])
        )
        idxs_of_one_bnds_frm_atm2 = idxs_of_other_bnds_frm_atm2[0]

        idx_of_othr_bnds_to_use.append(
            (idxs_of_one_bnds_frm_atm1, idxs_of_one_bnds_frm_atm2)
        )
        bonds_already_used.add(idxs_of_one_bnds_frm_atm1)
        bonds_already_used.add(idxs_of_one_bnds_frm_atm2)

    # Get all possible double-bond assignments (cis/trans).
    num = len(idx_of_othr_bnds_to_use)
    results = []
    if num == 0:
        # There are no unspecified double bonds, so just keep existing.
        results.append(mol)
        return results
    elif num == 1:
        # Only one unspecified double bond.
        all_chiral_options = [[True], [False]]
    else:
        # Many unspecified double bonds.
        starting = [[True], [False]]
        all_chiral_options = [[True], [False]]
        for i in range(num - 1):
            all_chiral_options = list(itertools.product(all_chiral_options, starting))
            all_chiral_options = [list(itertools.chain(c[0], c[1])) for c in all_chiral_options]

    # Let the user know.
    Utils.log(
        "\t" + mol.smiles(True) + " has " + str(int(len(all_chiral_options) / 2)) +
        " double bonds with unspecified stereochemistry."
    )

    # Pick a few to pursue futher. The aims to prevent a potential
    # combinatorial explosion.
    all_chiral_options = Utils.random_sample(
        all_chiral_options, max_variants_per_compound, ""
    )

    # Go through and consider each of the retained combinations.
    for chiral_options in all_chiral_options:
        # Make a copy of the original RDKit molecule.
        a_rd_mol = copy.copy(mol.rdkit_mol)

        # Make the assignment for this chiral_options.
        for idxs, chiral in zip(idx_of_othr_bnds_to_use, chiral_options):
            i1, i2 = idxs

            a_rd_mol.GetBondWithIdx(i1).SetBondDir(
                Chem.BondDir.ENDUPRIGHT
            )

            if chiral == True:
                a_rd_mol.GetBondWithIdx(i2).SetBondDir(
                    Chem.BondDir.ENDUPRIGHT
                )
            elif chiral == False:
                a_rd_mol.GetBondWithIdx(i2).SetBondDir(
                    Chem.BondDir.ENDDOWNRIGHT
                )

        # Assign the StereoChemistry. Required to actually set it.
        a_rd_mol.ClearComputedProps()
        Chem.AssignStereochemistry(a_rd_mol, force=True)

        # Make a new MyMol.MyMol object with the specified config.
        new_mol = MyMol.MyMol(a_rd_mol)

        if new_mol.can_smi != False and new_mol.can_smi != None:
            # Sometimes you get an error if there's a bad structure otherwise.

            # Add the new molecule to the list of results, if it does not have
            # a bizzare substructure.
            if not new_mol.remove_bizarre_substruc():
                new_mol.contnr_idx = mol.contnr_idx
                new_mol.name = mol.name
                new_mol.genealogy = mol.genealogy[:]
                new_mol.genealogy.append(
                    new_mol.smiles(True) + " (cis-trans isomerization)"
                )
                results.append(new_mol)

    # Return the results.
    return results
