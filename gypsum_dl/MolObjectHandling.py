##### MolObjectHandling.py

# Disable the unnecessary RDKit warnings
from rdkit import Chem, RDLogger

from gypsum_dl import utils

RDLogger.DisableLog("rdApp.*")


def check_sanitization(mol):
    """
    Given a Chem.rdchem.Mol this script will sanitize the molecule.
    It will be done using a series of try/except statements so that if it fails it will return a None
    rather than causing the outer script to fail.

    Nitrogen Fixing step occurs here to correct for a common RDKit valence error in which Nitrogens with
        with 4 bonds have the wrong formal charge by setting it to -1.
        This can be a place to add additional correcting features for any discovered common sanitation failures.

    Handled here so there are no problems later.

    Inputs:
    :param Chem.rdchem.Mol mol: an rdkit molecule to be sanitized
    Returns:
    :returns: Chem.rdchem.Mol mol: A sanitized rdkit molecule or None if it failed.
    """
    if mol is None:
        return None

    # easiest nearly everything should get through
    try:
        sanitize_string = Chem.SanitizeMol(
            mol,
            sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
            catchErrors=True,
        )
    except Exception:
        return None

    if sanitize_string.name == "SANITIZE_NONE":
        return mol

    # try to fix the nitrogen (common problem that 4 bonded Nitrogens improperly
    # lose their + charges)

    # The nitrogen fix guesses at formal charges, and the guess is only
    # justified if the molecule sanitizes afterwards. Work on a copy so that a
    # molecule this function ends up rejecting is not left in the caller's
    # hands carrying charges that were invented here.
    charges_before = [atom.GetFormalCharge() for atom in mol.GetAtoms()]
    candidate = Nitrogen_charge_adjustment(Chem.Mol(mol))
    if candidate is None:
        return None

    try:
        sanitize_string = Chem.SanitizeMol(
            candidate,
            sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
            catchErrors=True,
        )
    except Exception:
        return None

    # If any form of sanitation still fails (ie. KEKULIZE) then return None.
    if sanitize_string.name != "SANITIZE_NONE":
        return None

    # The adjustment decides a protonation state on the user's behalf, so it
    # has to be visible. This runs inside MyMol.__init__, before any genealogy
    # entry exists, which leaves the log as the only place to say so. The
    # comparison keeps the warning honest: the second sanitization pass can
    # succeed on its own, without any charge having been changed.
    if [atom.GetFormalCharge() for atom in candidate.GetAtoms()] != charges_before:
        utils.log(
            "\tWARNING: Adjusted a nitrogen formal charge to sanitize "
            + Chem.MolToSmiles(candidate)
        )

    return candidate


def try_deprotanation(sanitized_mol):
    """
    Given an already sanitize Chem.rdchem.Mol object, we will try to deprotanate the mol of all non-explicit
    Hs. If it fails it will return a None rather than causing the outer script to fail.

    Inputs:
    :param Chem.rdchem.Mol mol: an rdkit molecule already sanitized.
    Returns:
    :returns: Chem.rdchem.Mol mol_sanitized: an rdkit molecule with H's removed and sanitized.
                                            it returns None if H's can't be added or if sanitation fails
    """
    try:
        mol = Chem.RemoveHs(sanitized_mol, sanitize=False)
    except Exception:
        return None

    return check_sanitization(mol)


def try_reprotanation(sanitized_deprotanated_mol):
    """
    Given an already sanitize and deprotanate Chem.rdchem.Mol object, we will try to reprotanate the mol with
    implicit Hs. If it fails it will return a None rather than causing the outer script to fail.

    Inputs:
    :param Chem.rdchem.Mol sanitized_deprotanated_mol: an rdkit molecule already sanitized and deprotanated.
    Returns:
    :returns: Chem.rdchem.Mol mol_sanitized: an rdkit molecule with H's added and sanitized.
                                            it returns None if H's can't be added or if sanitation fails
    """

    if sanitized_deprotanated_mol is None:
        return None

    try:
        mol = Chem.AddHs(sanitized_deprotanated_mol)
    except Exception:
        mol = None

    return check_sanitization(mol)


#


def Nitrogen_charge_adjustment(mol):
    """
    When importing ligands with sanitation turned off, one can successfully import
    import a SMILES in which a Nitrogen (N) can have 4 bonds, but no positive charge.
    Any 4-bonded N lacking a positive charge will fail a sanitiation check.
        -This could be an issue with importing improper SMILES, reactions, or crossing a nuetral nitrogen
            with a side chain which adds an extra bond, but doesn't add the extra positive charge.

    To correct for this, this function will find all N atoms with a summed bond count of 4
    (ie. 4 single bonds;2 double bonds; a single and a triple bond; two single and a double bond)
    and set the formal charge of those N's to +1.

    RDkit treats aromatic bonds as a bond count of 1.5. But we will not try to correct for
    Nitrogens labeled as Aromatic. As precaution, any N which is aromatic is skipped in this function.

    Note that the adjustment is made in place; pass a copy if the original must
    be preserved.

    Inputs:
    :param Chem.rdchem.Mol mol: any rdkit mol
    Returns:
    :returns: Chem.rdchem.Mol mol: the same rdkit mol with the N's adjusted
    """
    if mol is None:
        return None
    # makes sure its an rdkit obj
    try:
        atoms = mol.GetAtoms()
    except Exception:
        return None

    for atom in atoms:
        if atom.GetAtomicNum() == 7:
            bonds = [bond.GetBondTypeAsDouble() for bond in atom.GetBonds()]
            # If aromatic skip as we do not want assume the charge.
            if 1.5 in bonds:
                continue
            # GetBondTypeAsDouble prints out 1 for single, 2.0 for double,
            # 3.0 for triple, 1.5 for AROMATIC but if AROMATIC WE WILL SKIP THIS ATOM
            num_bond_sums = sum(bonds)

            # Check if the octet is filled
            if num_bond_sums == 4.0:
                atom.SetFormalCharge(+1)
    return mol
