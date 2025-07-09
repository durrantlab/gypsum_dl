"""
The proccess_output definition determines which formats are saved to the
disk (output).
"""

from gypsum_dl.steps.io.to_html import web_2d_output
from gypsum_dl.steps.io.to_pdb import convert_sdfs_to_PDBs
from gypsum_dl.steps.io.to_sdf import save_to_sdf


def proccess_output(contnrs, params):
    """Proccess the molecular models in preparation for writing them to the
    disk."""

    # Unpack some variables.
    separate_output_files = params["separate_output_files"]
    output_folder = params["output_folder"]

    if params["add_html_output"] == True:
        # Write to an HTML file.
        web_2d_output(contnrs, output_folder)

    # Write to an SDF file.
    save_to_sdf(contnrs, params, separate_output_files, output_folder)

    # Also write to PDB files, if requested.
    if params["add_pdb_output"] == True:
        logger.info("Making PDB output files")
        convert_sdfs_to_PDBs(contnrs, output_folder)
