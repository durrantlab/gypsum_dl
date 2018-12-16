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
    :param max_variants_per_compound: To control the combinatorial explosion,
       only this number of variants (molecules) will be advanced to the next
       step.
    :type max_variants_per_compound: int
    :return: [description]
    :rtype: [type]
    """

    # For this to work, you need to have explicit hydrogens in place.
    mol.rdkit_mol = Chem.AddHs(mol.rdkit_mol)

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

    # Get a list of all the single bonds that come of each double-bond atom.
    all_sngl_bnd_idxs = set([])
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

        # idx_all_atms_in_dbl_bond.add(atom1.GetIdx())
        # idx_all_atms_in_dbl_bond.add(atom2.GetIdx())

        # Suffice it to say, RDKit does not deal with cis-trans isomerization
        # in an intuitive way...
        idxs_of_other_bnds_frm_atm1 = [b.GetIdx() for b in atom1.GetBonds()]
        idxs_of_other_bnds_frm_atm1.remove(dbl_bnd_idx)

        idxs_of_other_bnds_frm_atm2 = [b.GetIdx() for b in atom2.GetBonds()]
        idxs_of_other_bnds_frm_atm2.remove(dbl_bnd_idx)

        all_sngl_bnd_idxs |= set(idxs_of_other_bnds_frm_atm1)
        all_sngl_bnd_idxs |= set(idxs_of_other_bnds_frm_atm2)

    # Now come up with all possible up/down combinations for those bonds.
    all_sngl_bnd_idxs = list(all_sngl_bnd_idxs)
    all_chiral_options = list(
        itertools.product(
            [True, False],
            repeat=len(all_sngl_bnd_idxs)
        )
    )

    # # Get a list of all the single-bond pairs that come off each double-bond
    # # atom.
    # sngl_bnd_pairs = []
    # idx_all_atms_neigbr_dbl_bond = set([])
    # idx_all_atms_in_dbl_bond = set([])
    # for dbl_bnd_idx in unasignd_dbl_bnd_idxs:
    #     bond = mol.rdkit_mol.GetBondWithIdx(dbl_bnd_idx)

    #     atom1 = bond.GetBeginAtom()
    #     atom1_bonds = atom1.GetBonds()
    #     if len(atom1_bonds) == 1:
    #         # The only bond is the one you already know about. So don't save.
    #         continue

    #     atom2 = bond.GetEndAtom()
    #     atom2_bonds = atom2.GetBonds()
    #     if len(atom2_bonds) == 1:
    #         # The only bond is the one you already know about. So don't save.
    #         continue

    #     idx_all_atms_in_dbl_bond.add(atom1.GetIdx())
    #     idx_all_atms_in_dbl_bond.add(atom2.GetIdx())

    #     # Suffice it to say, RDKit does not deal with cis-trans isomerization
    #     # in an intuitive way...
    #     idxs_of_other_bnds_frm_atm1 = [b.GetIdx() for b in atom1.GetBonds()]
    #     idxs_of_other_bnds_frm_atm1.remove(dbl_bnd_idx)

    #     idxs_of_other_bnds_frm_atm2 = [b.GetIdx() for b in atom2.GetBonds()]
    #     idxs_of_other_bnds_frm_atm2.remove(dbl_bnd_idx)

    #     sngl_bnd_pairs.append(idxs_of_other_bnds_frm_atm1)
    #     sngl_bnd_pairs.append(idxs_of_other_bnds_frm_atm2)

    #     idx_all_atms_neigbr_dbl_bond |= set(idxs_of_other_bnds_frm_atm1)
    #     idx_all_atms_neigbr_dbl_bond |= set(idxs_of_other_bnds_frm_atm2)

    # # Now make up and down assignments based on those pairs. A recursive
    # # algorithm.
    # def recurse(assignments_so_far, sngl_bnd_pairs_left):
    #     """Recursively assigns bond directions.

    #     :param assignments_so_far: The assignments that have been made so far
    #        (growing list).
    #     :type assignments_so_far: list
    #     :param sngl_bnd_pairs_left: The atom pairs that have not yet been
    #        assigned up or down (diminishing list).
    #     :type sngl_bnd_pairs_left: list
    #     """

    #     if len(sngl_bnd_pairs_left) == 0:
    #         # Make sure to pass by value, not refernece.
    #         assignments_so_far = copy.copy(assignments_so_far)

    #         all_chiral_options.append(assignments_so_far)
    #         # End of recursion...
    #         return

    #     bond1_idx = sngl_bnd_pairs_left[0][0]
    #     bond1 = mol.rdkit_mol.GetBondWithIdx(bond1_idx)
    #     bond1.GetBeginAtom().GetIdx()

    #     import pdb; pdb.set_trace()

    #     bond2_idx = sngl_bnd_pairs_left[0][1]
    #     sngl_bnd_pairs_left = sngl_bnd_pairs_left[1:]

    #     if assignments_so_far[bond1_idx] == "up":
    #         # Make sure to pass by value, not refernece.
    #         assignments_so_far = copy.copy(assignments_so_far)

    #         # First atom already assigned, so second must be down.
    #         assignments_so_far[bond2_idx] = "down"
    #         recurse(assignments_so_far, sngl_bnd_pairs_left)
    #     elif assignments_so_far[bond1_idx] == "down":
    #         # Make sure to pass by value, not refernece.
    #         assignments_so_far = copy.copy(assignments_so_far)

    #         # First atom already assigned, so second must be up.
    #         assignments_so_far[bond2_idx] = "up"
    #         recurse(assignments_so_far, sngl_bnd_pairs_left)
    #     else:
    #         # Make sure to pass by value, not refernece.
    #         assignments_so_far1 = copy.copy(assignments_so_far)
    #         assignments_so_far2 = copy.copy(assignments_so_far)

    #         # Could go either way.
    #         assignments_so_far1[bond1_idx] = "up"
    #         assignments_so_far1[bond2_idx] = "down"
    #         recurse(assignments_so_far1, sngl_bnd_pairs_left)

    #         assignments_so_far2[bond1_idx] = "down"
    #         assignments_so_far2[bond2_idx] = "up"
    #         recurse(assignments_so_far2, sngl_bnd_pairs_left)

    # all_chiral_options = []
    # recurse({idx: "" for idx in idx_all_atms_neigbr_dbl_bond}, sngl_bnd_pairs)
    # for a in all_chiral_options:
    #     print(a)

    # # Throw out any bond that has an atom that only participate in that one
    # # bond (terminal alkene).
    # sngl_bnds_off_dble_infs = []
    # # bonds_already_used = set([])
    # for dbl_bnd_idx in unasignd_dbl_bnd_idxs:
    #     bond = mol.rdkit_mol.GetBondWithIdx(dbl_bnd_idx)

    #     atom1 = bond.GetBeginAtom()
    #     atom1_bonds = atom1.GetBonds()
    #     if len(atom1_bonds) == 1:
    #         # The only bond is the one you already know about. So don't save.
    #         continue

    #     atom2 = bond.GetEndAtom()
    #     atom2_bonds = atom2.GetBonds()
    #     if len(atom2_bonds) == 1:
    #         # The only bond is the one you already know about. So don't save.
    #         continue

    #     # Suffice it to say, RDKit does not deal with cis-trans isomerization
    #     # in an intuitive way...
    #     idxs_of_other_bnds_frm_atm1 = [b.GetIdx() for b in atom1.GetBonds()]
    #     idxs_of_other_bnds_frm_atm1.remove(dbl_bnd_idx)
    #     # idxs_of_other_bnds_frm_atm1 = list(
    #     #     set(idxs_of_other_bnds_frm_atm1) - bonds_already_used
    #     # )
    #     # idxs_of_one_bnds_frm_atm1 = idxs_of_other_bnds_frm_atm1[0]

    #     idxs_of_other_bnds_frm_atm2 = [b.GetIdx() for b in atom2.GetBonds()]
    #     idxs_of_other_bnds_frm_atm2.remove(dbl_bnd_idx)
    #     # idxs_of_other_bnds_frm_atm2 = list(
    #     #     set(idxs_of_other_bnds_frm_atm2) - bonds_already_used
    #     #     - set([idxs_of_one_bnds_frm_atm1])
    #     # )
    #     # idxs_of_one_bnds_frm_atm2 = idxs_of_other_bnds_frm_atm2[0]

    #     # sngl_bnds_off_dble_infs.append(
    #     #     (idxs_of_one_bnds_frm_atm1, idxs_of_one_bnds_frm_atm2)
    #     # )
    #     # bonds_already_used.add(idxs_of_one_bnds_frm_atm1)
    #     # bonds_already_used.add(idxs_of_one_bnds_frm_atm2)

    #     sngl_bnds_off_dble_infs.append({
    #         "fixed": idxs_of_other_bnds_frm_atm1,
    #         "movable": idxs_of_other_bnds_frm_atm2
    #     })

    # # Need to make sure no bond is labeled both fixed and movable.
    # already_movable_bonds = set([])
    # for sngl_bnds_off_dble_inf in sngl_bnds_off_dble_infs:
    #     fixed_bonds = sngl_bnds_off_dble_inf["fixed"]
    #     movable_bonds = sngl_bnds_off_dble_inf["movable"]

    #     # Have any of the bonds marked fixed been previously marked as
    #     # movable? If so, switch what is fixed and movable for this double
    #     # bond.
    #     if len(already_movable_bonds) != len(already_movable_bonds - set(fixed_bonds)):
    #         sngl_bnds_off_dble_inf["fixed"], sngl_bnds_off_dble_inf["movable"] = \
    #             sngl_bnds_off_dble_inf["movable"], sngl_bnds_off_dble_inf["fixed"]

    #     already_movable_bonds |= set(movable_bonds)

    # # Get all possible double-bond assignments (cis/trans).
    # num = len(sngl_bnds_off_dble_infs)
    # results = []
    # if num == 0:
    #     # There are no unspecified double bonds, so just keep existing.
    #     results.append(mol)
    #     return results
    # elif num == 1:
    #     # Only one unspecified double bond.
    #     all_chiral_options = [[True], [False]]
    # else:
    #     # Many unspecified double bonds.
    #     starting = [[True], [False]]
    #     all_chiral_options = [[True], [False]]
    #     for i in range(num - 1):
    #         all_chiral_options = list(itertools.product(all_chiral_options, starting))
    #         all_chiral_options = [list(itertools.chain(c[0], c[1])) for c in all_chiral_options]

    # Let the user know.
    Utils.log(
        "\t" + mol.smiles(True) + " has " + str(int(len(all_chiral_options) / 2)) +
        " double bonds with unspecified stereochemistry."
    )

    # Pick a few to pursue futher. The aims to prevent a potential
    # combinatorial explosion.
    # all_chiral_options = Utils.random_sample(
    #     all_chiral_options, max_variants_per_compound, ""
    # )

    import rdkit
    for bond in mol.rdkit_mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if bond.GetBondType() == rdkit.Chem.rdchem.BondType.SINGLE:
            c = "-"
        else:
            c = "="
        print("Bond", bond.GetIdx(), ":", atom1.GetSymbol(), c, atom2.GetSymbol())

    # Go through and consider each of the retained combinations.
    smiles_to_consider = set([])
    for chiral_options in all_chiral_options:
        # Make a copy of the original RDKit molecule.
        a_rd_mol = copy.copy(mol.rdkit_mol)
        # a_rd_mol = Chem.MolFromSmiles(mol.smiles())

        for bond_idx, direc in zip(all_sngl_bnd_idxs, chiral_options):
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

        # Add to list of ones to consider
        smiles_to_consider.add(
            Chem.MolToSmiles(a_rd_mol, isomericSmiles=True, canonical=True)
        )

    # Remove ones that don't have "/" or "\". These are not real enumerated ones.
    smiles_to_consider = [s for s in smiles_to_consider if "/" in s or "\\" in s]

    # Get the maximum number of / + \ in any string.
    cnts = [s.count("/") + s.count("\\") for s in smiles_to_consider]
    max_cnts = max(cnts)

    # Only keep those with that same max count. The others have double bonds
    # that remain unspecified.
    smiles_to_consider = [s[0] for s in zip(smiles_to_consider, cnts) if s[1] == max_cnts]

    print(cnts)
    print(smiles_to_consider)
    import pdb; pdb.set_trace()

    # # Go through and consider each of the retained combinations.
    # for chiral_options in all_chiral_options:
    #     # Make a copy of the original RDKit molecule.
    #     a_rd_mol = copy.copy(mol.rdkit_mol)
    #     # a_rd_mol = Chem.MolFromSmiles(mol.smiles())

    #     for bond_idx in chiral_options:
    #         direc = chiral_options[bond_idx]

    #         # Always done with reference to the atom in the double bond.

    #         if direc == "up":
    #             a_rd_mol.GetBondWithIdx(bond_idx).SetBondDir(
    #                 Chem.BondDir.ENDUPRIGHT
    #             )
    #         else:
    #             a_rd_mol.GetBondWithIdx(bond_idx).SetBondDir(
    #                 Chem.BondDir.ENDDOWNRIGHT
    #             )

        # Make the assignment for this chiral_options.
        # chiral_options[0] = True
        # for inf, chiral in zip(sngl_bnds_off_dble_infs, chiral_options):
        #     # i1, i2 = idxs
        #     i1 = inf["fixed"][0]  # The other one doesn't matter. It's
        #                           # direction is determined by the direction
        #                           # of the first one anyway.
        #     i2 = inf["movable"][0]
        #     print(i1, i2)

        #     a_rd_mol.GetBondWithIdx(i1).SetBondDir(
        #         Chem.BondDir.ENDUPRIGHT
        #     )

        #     if chiral == True:
        #         a_rd_mol.GetBondWithIdx(i2).SetBondDir(
        #             Chem.BondDir.ENDUPRIGHT
        #         )
        #     elif chiral == False:
        #         a_rd_mol.GetBondWithIdx(i2).SetBondDir(
        #             Chem.BondDir.ENDDOWNRIGHT
        #         )


        # print(chiral_options)  # If first one is true, doesn't work.
        # mm = copy.copy(a_rd_mol)
        # mm = Chem.RemoveHs(mm)
        # print(Chem.MolToSmiles(mm, isomericSmiles=True))
        # import pdb; pdb.set_trace()

    #     # Make a new MyMol.MyMol object with the specified config.
    #     new_mol = MyMol.MyMol(a_rd_mol)

    #     if new_mol.can_smi != False and new_mol.can_smi != None:
    #         # Sometimes you get an error if there's a bad structure otherwise.

    #         # Add the new molecule to the list of results, if it does not have
    #         # a bizzare substructure.
    #         if not new_mol.remove_bizarre_substruc():
    #             new_mol.contnr_idx = mol.contnr_idx
    #             new_mol.name = mol.name
    #             new_mol.genealogy = mol.genealogy[:]
    #             new_mol.genealogy.append(
    #                 new_mol.smiles(True) + " (cis-trans isomerization)"
    #             )
    #             results.append(new_mol)

    # # Return the results.
    # return results
