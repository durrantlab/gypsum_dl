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

import copy
import sys

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem.Draw import IPythonConsole #Needed to show molecules
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions #O
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign

def get_info_two_pdbs(mol_1_str, mol_2_str):
    if type(mol_1_str) == str:
        mol_1 = Chem.MolFromPDBFile(mol_1_str)
    if type(mol_2_str) == str:
        mol_2 = Chem.MolFromPDBFile(mol_2_str)

    ref = Chem.MolFromSmiles(Chem.MolToSmiles(mol_1))
    mol_A = AllChem.AssignBondOrdersFromTemplate(ref, mol_1)
    mol_B = AllChem.AssignBondOrdersFromTemplate(ref, mol_2)

    best_rmsd, num_rot_bond = get_RMSD_and_rot_bonds(mol_A, mol_B)

    return best_rmsd, num_rot_bond

def get_RMSD_and_rot_bonds(mol_ref, mol_test):
    # Align them
    best_rmsd = AllChem.GetBestRMS(mol_ref, mol_test)
    num_rot_bond=  rdkit.Chem.rdMolDescriptors.CalcNumRotatableBonds(mol_ref)

    return best_rmsd, num_rot_bond

if __name__ == "__main__":


    pdb_1 = str(sys.argv[1])
    pdb_2 = str(sys.argv[2])

    print(get_info_two_pdbs(pdb_1,pdb_2))
