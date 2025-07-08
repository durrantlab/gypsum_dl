"""
The proccess_output definition determines which formats are saved to the
disk (output).
"""

import __future__

from gypsum_dl import utils
from gypsum_dl.steps.io.SaveToPDB import convert_sdfs_to_PDBs
from gypsum_dl.steps.io.SaveToSDF import save_to_sdf
from gypsum_dl.steps.io.Web2DOutput import web_2d_output


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
        utils.log("\nMaking PDB output files\n")
        convert_sdfs_to_PDBs(contnrs, output_folder)
