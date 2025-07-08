"""
Saves output files to SDF.
"""

import os

from gypsum_dl import utils

try:
    from rdkit import Chem
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")


def save_to_sdf(contnrs, params, separate_output_files, output_folder):
    """Saves the 3D models to the disk as an SDF file.

    :param contnrs: A list of containers (container.MoleculeContainer).
    :type contnrs: list
    :param params: The parameters.
    :type params: dict
    :param separate_output_files: Whether save each molecule to a different
       file.
    :type separate_output_files: bool
    :param output_folder: The output folder.
    :type output_folder: str
    """

    # Save an empty molecule with the parameters.
    if separate_output_files == False:
        w = Chem.SDWriter(output_folder + os.sep + "gypsum_dl_success.sdf")
    else:
        w = Chem.SDWriter(output_folder + os.sep + "gypsum_dl_params.sdf")

    m = Chem.Mol()
    m.SetProp("_Name", "EMPTY MOLECULE DESCRIBING GYPSUM-DL PARAMETERS")
    for param in params:
        m.SetProp(param, str(params[param]))
    w.write(m)

    if separate_output_files == True:
        w.flush()
        w.close()

    # Also save the file or files containing the output molecules.
    utils.log("Saving molecules associated with...")
    for i, contnr in enumerate(contnrs):
        # Add the container properties to the rdkit_mol object so they get
        # written to the SDF file.
        contnr.add_container_properties()

        # Let the user know which molecule you're on.
        utils.log("\t" + contnr.orig_smi)

        # Save the file(s).
        if separate_output_files == True:
            # sdf_file = "{}{}__{}.pdb".format(output_folder + os.sep, slug(name), conformer_counter)
            sdf_file = f"{output_folder + os.sep}{utils.slug(contnr.name)}__input{contnr.contnr_idx_orig + 1}.sdf"
            w = Chem.SDWriter(sdf_file)
            # w = Chem.SDWriter(output_folder + os.sep + "output." + str(i + 1) + ".sdf")

        for m in contnr.mols:
            m.load_conformers_into_rdkit_mol()
            w.write(m.rdkit_mol)

        if separate_output_files == True:
            w.flush()
            w.close()

    if separate_output_files == False:
        w.flush()
        w.close()
