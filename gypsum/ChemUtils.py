"""The module includes definitions to manipulate the molecules."""

import __future__

import gypsum.Utils as Utils

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

def pick_lowest_enrgy_mols(mol_lst, num, thoroughness):
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
    :return: Returns a list of MyMol.MyMol, the best ones.
    :rtype: list
    """

    # Remove identical entries.
    mol_lst = list(set(mol_lst))

    # If the length of the mol_lst is less than num, just return them all.
    if len(mol_lst) <= num:
        return mol_lst

    # First, generate 3D structures. How many? num * thoroughness. mols_3d is
    # a list of gypsum MyMol.MyMol objects.
    mols_3d = Utils.random_sample(mol_lst, num * thoroughness, "")

    # Now get the energies
    data = []
    for i, mol in enumerate(mols_3d):
        mol.make_first_3d_conf_no_min()  # Make sure at least one conformer
                                         # exists.
        if len(mol.conformers) > 0:
            energy = mol.conformers[0].energy
            data.append((energy, i))

    data.sort()

    # Now keep only best top few.
    data = data[:num]

    # Keep just the mols there.
    new_mols_list = [mol_lst[d[1]] for d in data]

    # Return those molecules.
    return new_mols_list

def bst_for_each_contnr_no_opt(contnrs, mol_lst,
                               max_variants_per_compound, thoroughness,
                               crry_ovr_frm_lst_step_if_no_fnd=True):
    """Keep only the top few compound variants in each container, to prevent a
       combinatorial explosion. This is run periodically on the growing
       containers to keep them in check.

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
    """

    # Remove duplicate ligands from each container.
    for mol_cont in contnrs:
        mol_cont.remove_identical_mols_from_contnr()

    # Group the smiles by contnr_idx.
    data = Utils.group_mols_by_container_index(mol_lst)

    # Go through each container.
    for contnr_idx, contnr in enumerate(contnrs):
        none_generated = False

        # Pick just the lowest-energy conformers from the new candidates.
        # Possible a compound was eliminated early on, so doesn't exist.
        if contnr_idx in list(data.keys()):
            mols = data[contnr_idx]

            # Pick the lowest-energy molecules. Note that this creates a
            # conformation if necessary, but it is not minimized and so is not
            # computationally expensive.
            mols = pick_lowest_enrgy_mols(
                mols, max_variants_per_compound, thoroughness
            )

            if len(mols) > 0:
                # Now remove all previously determined mols for this
                # container.
                contnr.mols = []

                # Add in the lowest-energy conformers back to the smiles.
                for mol in mols:
                    contnr.add_mol(mol)
            else:
                none_generated = True
        else:
            none_generated = True

        # Now low-energy conformers were generated.
        if none_generated:
            if crry_ovr_frm_lst_step_if_no_fnd:
                # Just use previous ones.
                Utils.log(
                    "\tWARNING: Unable to find low-energy conformations: " +
                    contnr.orig_smi_deslt + " (" +
                    contnr.name + "). Keeping original " +
                    "conformers."
                )
            else:
                # Discard the conformation.
                Utils.log(
                    "\tWARNING: Unable to find low-energy conformations: " +
                    contnr.orig_smi_deslt + " (" +
                    contnr.name + "). Discarding conformer."
                )
                contnr.mols = []
