##### MolObjectHandling.py
import __future__

import rdkit
from rdkit import Chem
#Disable the unnecessary RDKit warnings
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

def check_sanitization(mol):
    """Given a rdkit.Chem.rdchem.Mol, this script will sanitize the molecule
    using a series of try/except statements so that if it fails, it will return
    a None rather than causing the outer script to fail.

    Nitrogen fixing occurs here also to correct for a common RDKit valence
    error in which nitrogens with with four bonds have the wrong formal charge.
    This function is a good place to add additional correcting features for any
    common sanitation failures discovered in the future.

    Handled here so there are no problems later.

    :param mol: An rdkit molecule to be sanitized
    :type mol: rdkit.Chem.rdchem.Mol
    :return: A sanitized rdkit molecule or None if it failed.
    :rtype: rdkit.Chem.rdchem.Mol | None
    """

    if mol is None:
        return None

    # Easiest nearly everything should get through.
    try:
        sanitize_string =  Chem.SanitizeMol(mol, sanitizeOps = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors = True)
    except:
        return None

    if sanitize_string.name == "SANITIZE_NONE":
        return mol
    else:
        # Try to fix the nitrogen (common problem that 4 bonded Nitrogens
        # improperly lose their + charges).
        mol = Nitrogen_charge_adjustment(mol)
        Chem.SanitizeMol(mol, sanitizeOps = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors = True)
        sanitize_string =  Chem.SanitizeMol(mol, sanitizeOps = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors = True)
        if sanitize_string.name == "SANITIZE_NONE":
            return mol

    # Run a sanitation Filter one more time in case something slipped through.
    # If there are any forms of sanitation which fail, i.e. KEKULIZE, then
    # return None.
    sanitize_string =  Chem.SanitizeMol(mol, sanitizeOps = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors = True)
    if sanitize_string.name != "SANITIZE_NONE":
        return None
    else:
        return mol

def try_deprotanation(sanitized_mol):
    """Given an already sanitized rdkit.Chem.rdchem.Mol object, we will try to
    deprotanate the mol of all non-explicit Hs. If it fails it will return a
    None rather than causing the outer script to fail.

    :param sanitized_mol: An rdkit molecule already sanitized.
    :type sanitized_mol: rdkit.Chem.rdchem.Mol
    :return: An rdkit molecule with H's removed and sanitized. It returns None
       if H's can't be added or if sanitation fails.
    :rtype: rdkit.Chem.rdchem.Mol | None
    """

    try:
        mol = Chem.RemoveHs(sanitized_mol, sanitize = False)
    except:
        return None

    mol_sanitized = check_sanitization(mol)

    return mol_sanitized

def try_reprotanation(sanitized_deprotanated_mol):
    """Given an already sanitized and deprotanated rdkit.Chem.rdchem.Mol
    object, we will try to reprotanate the mol with implicit Hs. If it fails,
    it will return a None rather than causing the outer script to fail.

    :param sanitized_deprotanated_mol: An rdkit molecule already sanitized and
       deprotanated.
    :type sanitized_deprotanated_mol: rdkit.Chem.rdchem.Mol
    :return: An rdkit molecule with H's added and sanitized. It returns None
       if H's can't be added or if sanitation fails.
    :rtype: rdkit.Chem.rdchem.Mol | None
    """

    if sanitized_deprotanated_mol is not None:
        try:
            mol = Chem.AddHs(sanitized_deprotanated_mol)
        except:
            mol = None

        mol_sanitized = check_sanitization(mol)
        return mol_sanitized
    else:
        return None
