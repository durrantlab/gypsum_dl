"""
Saves output files to SDF.
"""

import __future__

import contextlib
import os

from gypsum_dl import utils

try:
    from rdkit import Chem
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")

# Parallelizer stringifies to an address that changes every run, which makes
# byte-identical reruns impossible.
_PARAMS_NOT_WRITTEN = frozenset({"Parallelizer"})


def save_to_sdf(contnrs, params, separate_output_files, output_folder):
    """Saves the 3D models to the disk as an SDF file.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param params: The parameters.
    :type params: dict
    :param separate_output_files: Whether save each molecule to a different
       file.
    :type separate_output_files: bool
    :param output_folder: The output folder.
    :type output_folder: str
    """

    # The writer is tracked across the whole function so that an exception part
    # way through the molecules still leaves a flushed, closed file behind
    # rather than a truncated one. It is set back to None after each deliberate
    # close, so the handler at the bottom only ever finishes a writer that the
    # normal path did not reach.
    w = None

    try:
        # Save an empty molecule with the parameters.
        if separate_output_files == False:
            w = Chem.SDWriter(output_folder + os.sep + "gypsum_dl_success.sdf")
        else:
            # In separate-file mode (which MPI mode forces), every task would
            # otherwise open the same "gypsum_dl_params.sdf" for writing. Under
            # MPI all ranks run this concurrently against one path, interleaving
            # or truncating the record. One task holds exactly one container, so
            # qualify the params filename with that container's original index
            # to give each rank its own file.
            params_basename = "gypsum_dl_params"
            if contnrs:
                params_basename += f"__input{contnrs[0].contnr_idx_orig + 1}"
            w = Chem.SDWriter(f"{output_folder}{os.sep}{params_basename}.sdf")

        m = Chem.Mol()
        m.SetProp("_Name", "EMPTY MOLECULE DESCRIBING GYPSUM-DL PARAMETERS")
        for param in params:
            if param in _PARAMS_NOT_WRITTEN:
                continue
            m.SetProp(param, str(params[param]))
        w.write(m)

        if separate_output_files == True:
            w.flush()
            w.close()
            w = None

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
                if m.rdkit_mol is None:
                    # A molecule that failed sanitization or 3D optimization has
                    # nothing to write, and SDWriter.write(None) raises, which
                    # would take the rest of the output with it.
                    continue
                m.load_conformers_into_rdkit_mol()
                w.write(m.rdkit_mol)

            if separate_output_files == True:
                w.flush()
                w.close()
                w = None

        if separate_output_files == False:
            w.flush()
            w.close()
            w = None
    finally:
        if w is not None:
            # An exception is already on its way out of this function. Finish
            # the file so its records are not lost, but do not let a failure
            # here replace the error that caused it.
            with contextlib.suppress(Exception):
                w.flush()
                w.close()
