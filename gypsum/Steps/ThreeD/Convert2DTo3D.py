import __future__

import copy

import gypsum.parallelizer as parallelizer

import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

def mk3d(mol):
    show_error_msg = False

    if mol.rdkit_mol is None:
        show_error_msg = True
    else:
        if mol.remove_bizarre_substruc() == False:
            mol.makeMol3D()
            if len(mol.conformers) > 0:
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

def convert_2d_to_3d(contnrs, max_variants_per_compound, thoroughness, num_processors, multithread_mode, Parallelizer_obj):
    """
    Converts the 1D smiles strings into 3D small-molecule models.
    """

    Utils.log("Converting all molecules to 3D structures.")

    params = []
    for contnr in contnrs:
        for mol in contnr.mols:
            params.append([mol])

    tmp = Parallelizer_obj.run(params, mk3d, num_processors, multithread_mode)
    clear = parallelizer.strip_none(tmp)
    #clear = tmp
    ChemUtils.bst_for_each_contnr_no_opt(contnrs, clear, max_variants_per_compound, thoroughness, False)
