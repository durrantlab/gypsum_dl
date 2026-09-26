"""Unit tests for the RDKit sanitization and fragment helpers.

These helpers are the failure-handling layer for every malformed molecule
that enters the pipeline, so they are tested directly rather than only
through the end-to-end sample run.
"""

import pytest
from rdkit import Chem

import gypsum_dl.MolObjectHandling as MOH


def test_check_sanitization_none() -> None:
    assert MOH.check_sanitization(None) is None


def test_check_sanitization_handles_none_from_nitrogen_adjustment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression (B12): when Nitrogen_charge_adjustment returns None, the
    # subsequent SanitizeMol calls (outside any try) used to receive None and
    # raise. check_sanitization must instead return None cleanly. A neutral
    # quaternary nitrogen fails the first sanitize pass, so execution reaches
    # the adjustment branch.
    mol = Chem.MolFromSmiles("C[N](C)(C)C", sanitize=False)
    monkeypatch.setattr(MOH, "Nitrogen_charge_adjustment", lambda m: None)
    assert MOH.check_sanitization(mol) is None


def test_check_sanitization_accepts_valid_mol() -> None:
    mol = Chem.MolFromSmiles("c1ccccc1", sanitize=False)
    assert MOH.check_sanitization(mol) is not None


def test_check_sanitization_rejects_impossible_valence() -> None:
    mol = Chem.MolFromSmiles("C(C)(C)(C)(C)(C)C", sanitize=False)
    assert MOH.check_sanitization(mol) is None


def test_check_sanitization_charges_nitrogen_on_returned_mol() -> None:
    # A neutral quaternary nitrogen fails the first sanitize pass, so the
    # nitrogen fix runs. Whatever it corrects has to be present on the
    # molecule the caller gets back.
    fixed = MOH.check_sanitization(Chem.MolFromSmiles("C[N](C)(C)C", sanitize=False))
    assert fixed is not None
    charges = [a.GetFormalCharge() for a in fixed.GetAtoms() if a.GetAtomicNum() == 7]
    assert charges == [1]


def test_check_sanitization_does_not_charge_a_rejected_mol() -> None:
    # A neutral quaternary nitrogen alongside a hexavalent carbon: the nitrogen
    # fix applies but the molecule still fails to sanitize, so it is rejected.
    # A rejected molecule must not be handed back to the caller carrying the
    # +1 charge that the fix guessed at.
    mol = Chem.MolFromSmiles("C[N](C)(C)C.C(C)(C)(C)(C)C", sanitize=False)
    charges_before = [a.GetFormalCharge() for a in mol.GetAtoms()]
    assert all(c == 0 for c in charges_before)

    assert MOH.check_sanitization(mol) is None
    assert [a.GetFormalCharge() for a in mol.GetAtoms()] == charges_before


def test_check_sanitization_sanitizes_at_most_twice(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The rejection path used to call SanitizeMol four times: once up front,
    # twice in a row after the nitrogen fix (the first of the pair discarding
    # its result), and once more repeating the same call verbatim. Only the
    # first pass and one pass after the fix can change the outcome, and
    # sanitization runs on every molecule in the pipeline.
    calls: list[int] = []
    real_sanitize = Chem.SanitizeMol

    def counting_sanitize(*args: object, **kwargs: object) -> object:
        calls.append(1)
        return real_sanitize(*args, **kwargs)

    monkeypatch.setattr(MOH.Chem, "SanitizeMol", counting_sanitize)

    mol = Chem.MolFromSmiles("C(C)(C)(C)(C)(C)C", sanitize=False)
    assert MOH.check_sanitization(mol) is None
    assert len(calls) == 2


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


def test_remove_atoms_does_not_mutate_caller_list() -> None:
    # Regression (B13): remove_atoms aliased the caller's list and sorted it in
    # place, reordering the caller's data as a side effect. It must leave the
    # input list untouched.
    idx = [0, 2]
    MOH.remove_atoms(Chem.MolFromSmiles("CCO"), idx)
    assert idx == [0, 2]


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
