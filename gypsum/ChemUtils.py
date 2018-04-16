import sys
import gypsum.Utils as Utils

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)



def pick_best_mols(mol_lst, num, thoroughness):
    """
    Pick molecules with low energies.

    :param [MyMol.MyMol] mol_lst: The list of MyMol.MyMol objects.

    :param int num: The number of the lowest-energy ones to keep.

    :param int thoroughness: How many extra conformers to generate to pick num
               from. It's a factor.

    :returns: A list of MyMol.MyMol objects, the lowest energy ones generated.
    :rtype: :class:`str` ???
    """
    
    # Remove identical entries
    mol_lst = list(set(mol_lst))

    # If the length of the mol_lst is less than num, just return them all.
    if len(mol_lst) <= num:
        return mol_lst

    # First, generate 3D structures.
    # How many? num * thoroughness
    # mols_3d is a list of gypsum My.Mol instances #JAKE@
    mols_3d = Utils.random_sample(mol_lst, num * thoroughness, "")
    
    # Now get the energies
    data = []
    for i, mol in enumerate(mols_3d):
        mol.makeMol3D()  # make sure at least one conformer exists.
        if len(mol.conformers) > 0:
            energy = mol.conformers[0].energy
            data.append((energy, i))
    
    data.sort()

    # Now keep only best top few.
    data = data[:num]

    # Keep just the mols there
    new_mols_list = [mol_lst[d[1]] for d in data]

    # return those smiles
    return new_mols_list

def bst_for_each_contnr_no_opt(contnrs, mol_lst,
                                max_variants_per_compound, thoroughness,
                               crry_ovr_frm_lst_step_if_no_fnd=True):
    # smiles_tuple_data is (contnr_idx, MyMols)

    # Remove duplicate ligands from cluster.
    for mol_cont in contnrs:
        mol_cont.remove_identical_mols_from_container()

    # Group the smiles by contnr_idx
    data = Utils.group_mols_by_container_index(mol_lst)
    # Go through each contnr
    for contnr_idx, contnr in enumerate(contnrs):
        none_generated = False

        # Pick just the lowest-energy conformers from the new candidates.
        # Possible a compound was eliminated early on, so doesn't exist.
        if contnr_idx in data.keys():
            mols = data[contnr_idx]
            mols = pick_best_mols(
                mols, max_variants_per_compound, thoroughness
            )

            if len(mols) > 0:
                # Now remove all previously determined mols for this container
                contnr.mols = []

                # Add in the lowest-energy conformers back to the smiles.
                for mol in mols:
                    contnr.add_mol(mol)
            else:
                none_generated = True
        else:
            none_generated = True
        
        if none_generated:
            if crry_ovr_frm_lst_step_if_no_fnd:
                Utils.log(
                    "\tWARNING: Unable to find low-energy conformations: " +
                    contnr.orig_smi_deslt + " (" +
                    contnr.name + "). Keeping original " +
                    "conformers."
                )
            else:
                Utils.log(
                    "\tWARNING: Unable to find low-energy conformations: " +
                    contnr.orig_smi_deslt + " (" +
                    contnr.name + "). Discarding conformer."
                )
                contnr.mols = []
    