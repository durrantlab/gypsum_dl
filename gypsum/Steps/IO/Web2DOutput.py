# Copyright 2018 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Saves the output to an HTML file (2D images only). This is mostly for
debugging.
"""

import webbrowser
import os
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils

try:
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D
    from rdkit.Chem.Draw import PrepareMolForDrawing
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

def web_2d_output(contnrs, output_file):
    """Saves pictures of the models to an HTML file on disk. It can be viewed in
    a browser. Then opens a browser automatically to view them. This is mostly
    for debugging."""

    Utils.log("Saving html image of molecules associated with...")

    # Let's not parallelize it for now. This will rarely be used.
    f = open(output_file, 'w')
    for contnr in contnrs:
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

    # Open the browser to show the file.
    webbrowser.open("file://" + os.path.abspath(output_file))
