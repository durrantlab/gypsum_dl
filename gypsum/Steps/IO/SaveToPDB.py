import __future__

import glob
import sys
import os
from os.path import basename

import rdkit
import rdkit.Chem as Chem
#Disable the unnecessary RDKit warnings
rdkit.RDLogger.DisableLog('rdApp.*')

sys.path.append(os.path.join(os.path.abspath(os.path.dirname(__file__)),'gypsum'))
# from gypsum.gypsum import MolObjectHandling as MOH
 

def convert_sdfs_to_PDBs(contnrs, output_folder):
    """
    This will convert every conformer into a PDB file saved in the output_folder. 
    These .pdb files will be named after using the ID provided by the source file (or a modified name in the case of untitled or duplicate names)
        plus "__{}".format(i) where i is the conformer number
        
    Input:
    :param str gen_folder_path: Path of the folder for the current generation
    :param str SDFs_folder_path: Path of the folder with all of the 3D .sdf files to convert
    """
    for ligand in contnrs:
        ligand.add_container_properties()
        
        name = ligand.name
        mols = ligand.mols
        
        conformer_counter = 0
        for m in mols:
            pdb_file = "{}{}__{}.pdb".format(output_folder, name, conformer_counter)
            m.load_conformations_into_mol_3d()
            mol = m.rdkit_mol
            if mol == None:
                continue
            else:
                # Write to a PDB file
                Chem.MolToPDBFile(mol,pdb_file,flavor=4)
            conformer_counter += 1

#
            
            
