import __future__

from gypsum.Steps.IO.SaveToSDF import save_to_sdf
from gypsum.Steps.IO.SaveToPDB import convert_sdfs_to_PDBs
from gypsum.Steps.IO.Web2DOutput import web_2d_output

def proccess_output(contnrs, params):
    """
    Proccess outputing molecular models.
    """
    separate_output_files = params["separate_output_files"]
    output_file = params["output_file"]
    output_folder = params["output_folder"]

    if output_file.lower().endswith(".html"):
        web_2d_output(contnrs, output_file)
    else:
        save_to_sdf(contnrs, params, separate_output_files, output_file)


    if params["output_pdb"] == True:
        print("\nMaking PDB output files\n")
        convert_sdfs_to_PDBs(contnrs, output_folder)
