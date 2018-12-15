"""Saves output files to SDF."""

import __future__

import gypsum.Utils as Utils

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

def save_to_sdf(contnrs, params, separate_output_files, output_file):
    """Saves the 3D models to the disk as an SDF file.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param params: The parameters.
    :type params: dict
    :param separate_output_files: Whether save each molecule to a different
       file.
    :type separate_output_files: bool
    :param output_file: The output file (or partial file name if saving to
       multiple files).
    :type output_file: str
    """

    # Save an empty molecule with the parameters.
    if separate_output_files == False:
        w = Chem.SDWriter(output_file)
    else:
        w = Chem.SDWriter(output_file + ".params.sdf")

    m = Chem.Mol()
    m.SetProp("_Name", "EMPTY MOLECULE DESCRIBING GYPSUM PARAMETERS")
    for param in params:
        m.SetProp(param, str(params[param]))
    w.write(m)

    if separate_output_files == True:
        w.flush()
        w.close()

    # Also save the file or files containing the output molecules.
    Utils.log("Saving molecules associated with...")
    for i, contnr in enumerate(contnrs):
        # Add the container properties to the rdkit_mol object so they get
        # written to the SDF file.
        contnr.add_container_properties()

        # Let the user know which molecule you're on.
        Utils.log("\t" + contnr.orig_smi)

        # Save the file(s).
        if separate_output_files == True:
            w = Chem.SDWriter(output_file + "." + str(i + 1) + ".sdf")

        for m in contnr.mols:
            m.load_conformers_into_rdkit_mol()
            w.write(m.rdkit_mol)

        if separate_output_files == True:
            w.flush()
            w.close()

    if separate_output_files == False:
        w.flush()
        w.close()
