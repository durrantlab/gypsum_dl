"""
Saves the output to an HTML file (2D images only). This is mostly for
debugging.
"""

# import webbrowser
import os

from gypsum_dl import chem_utils, utils

try:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import PrepareMolForDrawing, rdMolDraw2D
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")


def web_2d_output(contnrs, output_folder):
    """Saves pictures of the models to an HTML file on disk. It can be viewed in
    a browser. Then opens a browser automatically to view them. This is mostly
    for debugging."""

    utils.log("Saving html image of molecules associated with...")

    # Let's not parallelize it for now. This will rarely be used.
    html_file = output_folder + os.sep + "gypsum_dl_success.html"
    with open(html_file, "w") as f:
        for contnr in contnrs:
            utils.log("\t" + contnr.orig_smi)
            for mol in contnr.mols:
                # See
                # http://rdkit.org/docs/source/rdkit.Chem.rdmolops.html#rdkit.Chem.rdmolops.RemoveHs
                # I think in older versions of rdkit (e.g., 2016.09.2), RemoveHs
                # would remove hydrogens, even if that make double bonds
                # ambiguous. Not so in newer versions (e.g., 2018.03.4). So if
                # your double-bonded nitrogen doesn't have its hydrogen attached,
                # and you're using an older version of rdkit, don't worry about
                # it. The cis/trans info is still there.
                mol2 = Chem.RemoveHs(mol.rdkit_mol)
                # mol2 = mol.rdkit_mol

                mol2 = PrepareMolForDrawing(mol2, addChiralHs=True, wedgeBonds=True)
                rdDepictor.Compute2DCoords(mol2)
                drawer = rdMolDraw2D.MolDraw2DSVG(200, 200)
                drawer.DrawMolecule(mol2)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                f.write(
                    '<div style="float: left; width:200px; height: 220px;" title="'
                    + mol.name
                    + '">'
                    + '<div style="width: 200px; height: 200px;">'
                    + svg.replace("svg:", "")
                    + "</div>"
                    + '<div style="width: 200px; height: 20px;">'
                    + "<small><center>"
                    + mol.smiles(True)
                    + "</center></small>"
                    + "</div>"
                    + "</div>"
                )

    # Open the browser to show the file.
    # webbrowser.open("file://" + os.path.abspath(html_file))
