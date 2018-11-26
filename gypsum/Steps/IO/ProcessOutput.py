import __future__

from gypsum.Steps.IO.SaveToSDF import save_to_sdf
from gypsum.Steps.IO.Web2DOutput import web_2d_output

def proccess_output(contnrs, params):
    """
    Proccess outputing molecular models.
    """
    separate_output_files = params["separate_output_files"]
    output_file = params["output_file"]

    if output_file.lower().endswith(".html"):
        web_2d_output(contnrs, output_file)
    else:
        save_to_sdf(contnrs, params, separate_output_files, output_file)
