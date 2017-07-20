from gypsum.Steps.IO.SaveToSDF import save_to_sdf
from gypsum.Steps.IO.Web2DOutput import web_2d_output

def proccess_output(self):
    """
    Proccess outputing molecular models.
    """
    if self.params["output_file"].lower().endswith(".html"):
        web_2d_output(self)
    else:
        save_to_sdf(self)
