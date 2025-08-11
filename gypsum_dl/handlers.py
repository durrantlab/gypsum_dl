from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")  # type: ignore


def check_sanitization(mol: Chem.Mol) -> Chem.Mol | None:
    """
    Given a Chem.rdchem.Mol this script will sanitize the molecule.
    It will be done using a series of try/except statements so that if it fails it will return a None
    rather than causing the outer script to fail.

    Nitrogen Fixing step occurs here to correct for a common RDKit valence error in which Nitrogens with
    with 4 bonds have the wrong formal charge by setting it to -1.
    This can be a place to add additional correcting features for any discovered common sanitation failures.

    Handled here so there are no problems later.

    Args:
        mol: an rdkit molecule to be sanitized

    Returns:
        A sanitized rdkit molecule or None if it failed.
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
    mol = Nitrogen_charge_adjustment(mol)
    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
        catchErrors=True,
    )
    sanitize_string = Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
        catchErrors=True,
    )
    if sanitize_string.name == "SANITIZE_NONE":
        return mol

    # run a  sanitation Filter 1 more time incase something slipped through ie.
    # if there are any forms of sanition which fail ie. KEKULIZE then return
    # None
    sanitize_string = Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
        catchErrors=True,
    )

    return None if sanitize_string.name != "SANITIZE_NONE" else mol


def handleHs(mol: Chem.Mol, protanate_step: bool):
    """
    Given a Chem.rdchem.Mol this script will sanitize the molecule, remove all non-explicit H's
    and add back on all implicit H's. This is to control for any discrepencies in the smiles strings or presence and
    absense of H's.
    If it fails it will return a None rather than causing the outer script to fail. Handled here so there are no problems later.

    Inputs:
        sanitized_deprotanated_mol: an rdkit molecule already sanitized and deprotanated.
        protanate_step: True if mol needs to be protanated; False if deprotanated
            Note if Protanated, SmilesMerge takes up to 10times longer

    Returns:
        An rdkit molecule with H's handled (either added or removed) and sanitized.
            it returns None if H's can't be added or if sanitation fails
    """
    mol = check_sanitization(mol)
    if mol is None:
        # mol failed Sanitation
        return None

    mol = try_deprotanation(mol)
    if mol is None:
        # mol failed deprotanation
        return None

    if protanate_step is True:
        # PROTANTION IS ON
        mol = try_reprotanation(mol)
        if mol is None:
            # mol failed reprotanation
            return None

    return mol


def try_deprotanation(sanitized_mol: Chem.Mol) -> Chem.Mol | None:
    """
    Given an already sanitize Chem.rdchem.Mol object, we will try to deprotanate the mol of all non-explicit
    Hs. If it fails it will return a None rather than causing the outer script to fail.

    Args:
        mol: an rdkit molecule already sanitized.

    Returns:
        An rdkit molecule with H's removed and sanitized.
            it returns None if H's can't be added or if sanitation fails.
    """
    try:
        mol = Chem.RemoveHs(sanitized_mol, sanitize=False)
    except Exception:
        return None

    return check_sanitization(mol)


def try_reprotanation(sanitized_deprotanated_mol: Chem.Mol) -> Chem.Mol | None:
    """
    Given an already sanitize and deprotanate Chem.rdchem.Mol object, we will try to reprotanate the mol with
    implicit Hs. If it fails it will return a None rather than causing the outer script to fail.

    Args:
        sanitized_deprotanated_mol: an rdkit molecule already sanitized and deprotanated.

    Returns:
        An rdkit molecule with H's added and sanitized.
            it returns None if H's can't be added or if sanitation fails
    """

    if sanitized_deprotanated_mol is None:
        return None

    try:
        mol = Chem.AddHs(sanitized_deprotanated_mol)
    except Exception:
        mol = None

    return check_sanitization(mol)


def remove_atoms(mol: Chem.Mol, list_of_idx_to_remove: list[int]) -> Chem.Mol | None:
    """
    This function removes atoms from an rdkit mol based on
    a provided list. The RemoveAtom function in Rdkit requires
    converting the mol to an more editable version of the rdkit mol
    object (Chem.EditableMol).

    Args:
        mol: any rdkit mol
        list_of_idx_to_remove: a list of idx values to remove
            from mol
    Returns:
        the rdkit mol as input but with
            the atoms from the list removed
    """

    if mol is None:
        return None

    try:
        atomsToRemove = list_of_idx_to_remove
        atomsToRemove.sort(reverse=True)
    except Exception:
        return None

    try:
        em1 = Chem.EditableMol(mol)
        for atom in atomsToRemove:
            em1.RemoveAtom(atom)

        return em1.GetMol()
    except Exception:
        return None


def Nitrogen_charge_adjustment(mol: Chem.Mol) -> Chem.Mol | None:
    """
    When importing ligands with sanitation turned off, one can successfully import
    import a SMILES in which a Nitrogen (N) can have 4 bonds, but no positive charge.
    Any 4-bonded N lacking a positive charge will fail a sanitiation check.
    This could be an issue with importing improper SMILES, reactions, or crossing a neutral nitrogen
    with a side chain which adds an extra bond, but doesn't add the extra positive charge.

    To correct for this, this function will find all N atoms with a summed bond count of 4
    (ie. 4 single bonds;2 double bonds; a single and a triple bond; two single and a double bond)
    and set the formal charge of those N's to +1.

    RDkit treats aromatic bonds as a bond count of 1.5. But we will not try to correct for
    Nitrogens labeled as Aromatic. As precaution, any N which is aromatic is skipped in this function.

    Args:
        mol: any rdkit mol

    Returns:
        the same rdkit mol with the N's adjusted
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


def check_for_unassigned_atom(mol: Chem.Mol) -> Chem.Mol | None:
    """
    Check there isn't a missing atom group ie. '*'
    A '*' in a SMILES string is an atom with an atomic num of 0
    """
    if mol is None:
        return None

    try:
        atoms = mol.GetAtoms()
    except Exception:
        return None

    for atom in atoms:
        if atom.GetAtomicNum() == 0:
            return None
    return mol


def handle_frag_check(mol: Chem.Mol) -> Chem.Mol | None:
    """
    This will take a RDKit Mol object. It will check if it is fragmented.
    If it has fragments it will return the largest of the two fragments.
    If it has no fragments it will return the molecule on harmed.
    """
    if mol is None:
        return None

    try:
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    except Exception:
        return None

    if len(frags) == 1:
        return mol
    else:
        frag_info_list = []
        frag_index = 0
        for frag in frags:
            # Check for unassigned breaks ie. a '*'
            frag = check_for_unassigned_atom(frag)
            if frag is None:
                frag_index = frag_index + 1
                continue
            else:
                num_atoms = frag.GetNumAtoms()
                frag_info = [frag_index, num_atoms]
                frag_info_list.append(frag_info)
                frag_index = frag_index + 1
        if len(frag_info_list) == 0:
            return None
        # Get the largest Fragment
        frag_info_list.sort(key=lambda x: float(x[-1]), reverse=True)
        largest_frag_idx = frag_info_list[0][0]
        largest_frag = frags[largest_frag_idx]
        return largest_frag
