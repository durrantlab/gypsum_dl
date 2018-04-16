import copy
import sys

import gypsum.Multiprocess as mp
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)

def mk3d(mol):
    show_error_msg = False

    if mol.rdkit_mol is None:
        show_error_msg = True
    else:
        if mol.crzy_substruc() == False:
            mol.makeMol3D()
            if mol.GetNumConformers() > 0:
                mol.genealogy.append(
                    mol.smiles(True) + " (3D coordinates assigned)"
                )
                return mol
            else:
                show_error_msg = True

    if show_error_msg:
        Utils.log(
            "\tWarning: Could not generate 3D geometry for " +
            str(mol.smiles()) + " (" + mol.name + "). Molecule " +
            "discarded."
        )

def convert_2d_to_3d(contnrs, max_variants_per_compound, thoroughness, num_processors):
    """
    Converts the 1D smiles strings into 3D small-molecule models.
    """

    Utils.log("Converting all molecules to 3D structures.")

    params = []
    for contnr in contnrs:
        for mol in contnr.mols:
            params.append(mol)

    tmp = mp.MultiThreading(params, num_processors, mk3d)

    clear = mp.strip_none(tmp)
    #clear = tmp
    ChemUtils.bst_for_each_contnr_no_opt(contnrs, clear, max_variants_per_compound, thoroughness, False)
