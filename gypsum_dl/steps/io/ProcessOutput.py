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

    # Write to an SDF file. This is the primary deliverable, so it comes first:
    # the optional formats below depict every variant, and a single variant
    # that survives the pipeline but cannot be drawn must not cost the user the
    # output of an entire run.
    save_to_sdf(contnrs, params, separate_output_files, output_folder)

    # Also write to PDB files, if requested.
    if params["add_pdb_output"] == True:
        utils.log("\nMaking PDB output files\n")
        convert_sdfs_to_PDBs(contnrs, output_folder)

    # Write to an HTML file, if requested. This one is for debugging, so it is
    # not worth ending the run over.
    if params["add_html_output"] == True:
        try:
            web_2d_output(contnrs, output_folder)
        except Exception:
            utils.log(
                "WARNING: Could not write the HTML (2D depiction) output. The "
                "SDF output was written successfully."
            )
