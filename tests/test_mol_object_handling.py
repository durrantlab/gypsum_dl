"""Unit tests for the RDKit sanitization and fragment helpers.

These helpers are the failure-handling layer for every malformed molecule
that enters the pipeline, so they are tested directly rather than only
through the end-to-end sample run.
"""

from rdkit import Chem

import gypsum_dl.MolObjectHandling as MOH


def test_check_sanitization_none() -> None:
    assert MOH.check_sanitization(None) is None


def test_check_sanitization_accepts_valid_mol() -> None:
    mol = Chem.MolFromSmiles("c1ccccc1", sanitize=False)
    assert MOH.check_sanitization(mol) is not None


def test_check_sanitization_rejects_impossible_valence() -> None:
    mol = Chem.MolFromSmiles("C(C)(C)(C)(C)(C)C", sanitize=False)
    assert MOH.check_sanitization(mol) is None


def test_try_deprotanation_removes_explicit_hydrogens() -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    deprotonated = MOH.try_deprotanation(mol)
    assert deprotonated is not None
    assert deprotonated.GetNumAtoms() == 3


def test_try_reprotanation_adds_hydrogens() -> None:
    reprotonated = MOH.try_reprotanation(Chem.MolFromSmiles("CCO"))
    assert reprotonated is not None
    assert reprotonated.GetNumAtoms() == 9


def test_try_reprotanation_none() -> None:
    assert MOH.try_reprotanation(None) is None


def test_handle_hs_without_protonation_step() -> None:
    mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    result = MOH.handleHs(mol, False)
    assert result is not None
    assert result.GetNumAtoms() == 3


def test_handle_hs_with_protonation_step() -> None:
    result = MOH.handleHs(Chem.MolFromSmiles("CCO"), True)
    assert result is not None
    assert result.GetNumAtoms() == 9


def test_handle_hs_none() -> None:
    assert MOH.handleHs(None, True) is None


def test_remove_atoms_drops_requested_indices() -> None:
    trimmed = MOH.remove_atoms(Chem.MolFromSmiles("CCO"), [2])
    assert trimmed is not None
    assert trimmed.GetNumAtoms() == 2


def test_remove_atoms_none_mol() -> None:
    assert MOH.remove_atoms(None, [0]) is None


def test_remove_atoms_unsortable_index_list() -> None:
    assert MOH.remove_atoms(Chem.MolFromSmiles("CCO"), 0) is None


def test_nitrogen_charge_adjustment_charges_quaternary_nitrogen() -> None:
    mol = Chem.MolFromSmiles("C[N](C)(C)C", sanitize=False)
    adjusted = MOH.Nitrogen_charge_adjustment(mol)
    assert adjusted is not None
    charges = [a.GetFormalCharge() for a in adjusted.GetAtoms() if a.GetAtomicNum() == 7]
    assert charges == [1]


def test_nitrogen_charge_adjustment_skips_aromatic_nitrogen() -> None:
    adjusted = MOH.Nitrogen_charge_adjustment(Chem.MolFromSmiles("c1ccncc1"))
    assert adjusted is not None
    charges = [a.GetFormalCharge() for a in adjusted.GetAtoms() if a.GetAtomicNum() == 7]
    assert charges == [0]


def test_nitrogen_charge_adjustment_none() -> None:
    assert MOH.Nitrogen_charge_adjustment(None) is None


def test_check_for_unassigned_atom_rejects_dummy() -> None:
    assert MOH.check_for_unassigned_atom(Chem.MolFromSmiles("*CC")) is None


def test_check_for_unassigned_atom_accepts_real_atoms() -> None:
    mol = Chem.MolFromSmiles("CCO")
    assert MOH.check_for_unassigned_atom(mol) is mol


def test_check_for_unassigned_atom_none() -> None:
    assert MOH.check_for_unassigned_atom(None) is None


def test_handle_frag_check_single_fragment() -> None:
    mol = Chem.MolFromSmiles("CCO")
    assert MOH.handle_frag_check(mol) is mol


def test_handle_frag_check_keeps_largest_fragment() -> None:
    largest = MOH.handle_frag_check(Chem.MolFromSmiles("CCCCCCO.C"))
    assert largest is not None
    assert largest.GetNumAtoms() == 7


def test_handle_frag_check_skips_dummy_fragments() -> None:
    largest = MOH.handle_frag_check(Chem.MolFromSmiles("*.CCO"))
    assert largest is not None
    assert largest.GetNumAtoms() == 3


def test_handle_frag_check_all_fragments_unassigned() -> None:
    assert MOH.handle_frag_check(Chem.MolFromSmiles("*.*")) is None


def test_handle_frag_check_none() -> None:
    assert MOH.handle_frag_check(None) is None