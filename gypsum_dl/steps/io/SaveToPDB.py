"""
Contains the function for saving the output to PDB files.
"""

import os
import sys

import rdkit.Chem as Chem

from gypsum_dl import utils

sys.path.append(os.path.join(os.path.abspath(os.path.dirname(__file__)), "gypsum_dl"))


def convert_sdfs_to_PDBs(contnrs, output_folder):
    """This will convert every conformer into a PDB file, which is saved in
       the output_folder.

    :param contnrs: A list of containers (container.MoleculeContainer).
    :type contnrs: list
    :param output_folder: The name of the output folder.
    :type output_folder: str
    """

    # Go through each container.
    for contnr in contnrs:
        contnr.add_container_properties()

        # Get the molecule name and associated variants.
        name = contnr.name
        mols = contnr.mols

        # Got through the variants.
        for i, m in enumerate(mols):
            pdb_file = f"{output_folder + os.sep}{utils.slug(name)}__input{contnr.contnr_idx_orig + 1}__variant{i + 1}.pdb"

            # Get the conformers into the rdkit_mol object.
            m.load_conformers_into_rdkit_mol()
            mol = m.rdkit_mol
            if mol is None:
                continue

            # Write conformers to a PDB file.
            Chem.MolToPDBFile(mol, pdb_file, flavor=32)

            # Add header to PDB file with original SMILES and final SMILES
            printout = f"REMARK Original SMILES string: {m.orig_smi}\nREMARK Final SMILES string: {m.standardize_smiles()}\n"
            with open(pdb_file) as f:
                printout += f.read()
            with open(pdb_file, "w") as f:
                f.write(printout)
            printout = ""
