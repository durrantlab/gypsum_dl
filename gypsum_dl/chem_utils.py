"""The module includes definitions to manipulate the molecules."""

from typing import TYPE_CHECKING

from loguru import logger
from rdkit import Chem

from gypsum_dl import utils

if TYPE_CHECKING:
    from gypsum_dl.models import Molecule, MoleculeContainer


def pick_lowest_enrgy_mols(
    mol_lst: list["Molecule"], num: int, thoroughness: int
) -> list["Molecule"]:
    """Pick molecules with low energies. If necessary, the definition also
       makes a conformer without minimization (so not too computationally
       expensive).

    Args:
        mol_lst: The list of Molecule objects.
        num: The number of the lowest-energy ones to keep.
        thoroughness: How many molecules to generate per variant (molecule)
            retained, for evaluation. For example, perhaps you want to advance five
            molecules (max_variants_per_compound = 5). You could just generate five
            and advance them all. Or you could generate ten and advance the best
            five (so thoroughness = 2). Using thoroughness > 1 increases the
            computational expense, but it also increases the chances of finding good
            molecules.

    Returns:
        Returns a list of Molecule, the best ones.
    """

    # Remove identical entries.
    mol_lst = list(set(mol_lst))

    # If the length of the mol_lst is less than num, just return them all.
    if len(mol_lst) <= num:
        return mol_lst

    # First, generate 3D structures. How many? num * thoroughness. mols_3d is
    # a list of Gypsum-DL Molecule objects.
    mols_3d = utils.random_sample(mol_lst, num * thoroughness, "")

    # Now get the energies
    data = []
    for i, mol in enumerate(mols_3d):
        mol.make_first_3d_conf_no_min()  # Make sure at least one conformer exists.

        if len(mol.conformers) > 0:
            energy = mol.conformers[0].energy
            data.append((energy, i))

    data.sort()

    # Now keep only best top few.
    data = data[:num]

    # Keep just the mols there.
    return [mol_lst[d[1]] for d in data]


def remove_highly_charged_molecules(mol_lst: list["Molecule"]) -> list["Molecule"]:
    """Remove molecules that are highly charged.

    Args:
        mol_lst: The list of molecules to consider.

    Returns:
        A list of molecules that are not too charged.
    """

    # First, find the molecule that is closest to being neutral.
    charges = [Chem.GetFormalCharge(mol.rdkit_mol) for mol in mol_lst]
    abs_charges = [abs(c) for c in charges]
    idx_of_closest_to_neutral = abs_charges.index(min(abs_charges))
    charge_closest_to_neutral = charges[idx_of_closest_to_neutral]

    # Now create a new mol list, where the charges deviation from the most
    # neutral by no more than 4. Note that this used to be 2, but I increased
    # it to 4 to accommodate ATP.
    new_mol_lst = []
    for i, charge in enumerate(charges):
        if abs(charge - charge_closest_to_neutral) <= 4:
            new_mol_lst.append(mol_lst[i])
        else:
            logger.warning(
                "Discarding highly charged form: " + mol_lst[i].smiles() + "."
            )

    return new_mol_lst


def bst_for_each_contnr_no_opt(
    contnrs: list["MoleculeContainer"],
    mol_lst: list["Molecule"],
    max_variants_per_compound: int,
    thoroughness: int,
    crry_ovr_frm_lst_step_if_no_fnd: bool = True,
) -> None:
    """Keep only the top few compound variants in each container, to prevent a
    combinatorial explosion. This is run periodically on the growing
    containers to keep them in check.

    Args:
        contnrs: A list of containers (container.MoleculeContainer).
        mol_lst: The list of Molecule objects.
        max_variants_per_compound: To control the combinatorial explosion,
            only this number of variants (molecules) will be advanced to the next
            step.
        thoroughness: How many molecules to generate per variant (molecule)
            retained, for evaluation. For example, perhaps you want to advance five
            molecules (max_variants_per_compound = 5). You could just generate five
            and advance them all. Or you could generate ten and advance the best
            five (so thoroughness = 2). Using thoroughness > 1 increases the
            computational expense, but it also increases the chances of finding good
            molecules.
        crry_ovr_frm_lst_step_if_no_fnd: If it can't find any low-energy
            conformers, determines whether to just keep the old ones. Defaults to
            True.
    """

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

            # Pick the lowest-energy molecules. Note that this creates a
            # conformation if necessary, but it is not minimized and so is not
            # computationally expensive.
            mols = pick_lowest_enrgy_mols(mols, max_variants_per_compound, thoroughness)

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

        # No low-energy conformers were generated.
        if none_generated:
            if crry_ovr_frm_lst_step_if_no_fnd:
                # Just use previous ones.
                logger.warning(
                    "Unable to find low-energy conformations: "
                    + contnr.orig_smi_deslt
                    + " ("
                    + contnr.name
                    + "). Keeping original "
                    + "conformers."
                )
            else:
                # Discard the conformation.
                logger.warning(
                    "Unable to find low-energy conformations: "
                    + contnr.orig_smi_deslt
                    + " ("
                    + contnr.name
                    + "). Discarding conformer."
                )
                contnr.mols = []


def uniq_mols_in_list(mol_lst: list[Chem.Mol]) -> list[str]:
    # You need to make new molecules to get it to work.
    # new_smiles = [m.smiles() for m in self.mols]
    # new_mols = [Chem.MolFromSmiles(smi) for smi in new_smiles]
    # new_can_smiles = [Chem.MolToSmiles(new_mol, isomericSmiles=True, canonical=True) for new_mol in new_mols]

    can_smiles_already_set = set([])
    uniq_mols = []
    for m in mol_lst:
        smi: str = m.smiles()
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
