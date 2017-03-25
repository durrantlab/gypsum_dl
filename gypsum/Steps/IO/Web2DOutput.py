import webbrowser
import os
import sys
from ... import Utils
from ... import ChemUtils

try:
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem.Draw import PrepareMolForDrawing
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)
    
def web_2d_output(self):
    """
    Saves pictures of the models to an HTML file on disk. It can be viewed in
    a browser. Then opens a browser automatically to view them.
    """

    Utils.log("Saving html image of molecules associated with...")

    # Let's not parallelize it for now.
    f = open(self.params["output_file"], 'w')
    for contnr in self.contnrs:
        Utils.log("\t" + contnr.orig_smi)
        for mol in contnr.mols:
            mol2 = Chem.RemoveHs(mol.rdkit_mol)
            mol2 = PrepareMolForDrawing(mol2, addChiralHs=True, wedgeBonds=True)
            rdDepictor.Compute2DCoords(mol2)
            drawer = rdMolDraw2D.MolDraw2DSVG(200,200)
            drawer.DrawMolecule(mol2)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()
            f.write(
                '<div style="float: left; width:200px; height: 220px;">' + 
                    '<div style="width: 200px; height: 200px;">' + 
                        svg.replace("svg:", "") + 
                    '</div>' +
                    '<div style="width: 200px; height: 20px;">' + 
                        '<small><center>' + mol.smiles(True) + '</center></small>' +
                    '</div>' + 
                '</div>')
    f.close()

    webbrowser.open("file://" + os.path.abspath(self.params["output_file"]))
