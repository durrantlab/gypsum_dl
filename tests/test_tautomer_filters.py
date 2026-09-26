"""Tests for the tautomer rejection filters.

`tauts_no_change_hs_to_cs_unless_alpha_to_carbnyl` is not wired into
`prepare_smiles`, so it and its worker are only reachable from here.
"""

import pytest

from gypsum_dl import MyMol
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.steps.smiles import MakeTautomers
from gypsum_dl.steps.smiles.MakeTautomers import (
    parallel_check_carbon_hydrogens,
    parallel_check_chiral_centers,
    parallel_check_nonarom_rings,
    parallel_make_taut,
    tauts_no_break_arom_rngs,
    tauts_no_change_hs_to_cs_unless_alpha_to_carbnyl,
    tauts_no_elim_chiral,
)


def _taut(smiles: str, name: str) -> MyMol.MyMol:
    """Build a candidate tautomer tagged for container zero.

    Args:
        smiles: SMILES string for the tautomer.
        name: Ligand name used in the rejection log messages.

    Returns:
        The tagged MyMol instance.
    """
    mol = MyMol.MyMol(smiles)
    mol.name = name
    mol.contnr_idx = 0
    return mol


def test_parallel_make_taut_returns_none_when_rdkit_mol_unsanitizable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: parallel_make_taut called Chem.RemoveHs before its None
    # check, so an unsanitizable molecule (rdkit_mol is None) raised inside
    # RemoveHs instead of being dropped gracefully.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.mols.append(_taut("CCO", "ethanol"))

    class _NoneMol:
        rdkit_mol = None

    monkeypatch.setattr(MakeTautomers.MyMol, "MyMol", lambda *a, **k: _NoneMol())
    assert parallel_make_taut(contnr, 0, 1) is None


def test_parallel_check_nonarom_rings_keeps_matching_tautomer() -> None:
    contnr = MolContainer("C1CCCCC1", "cyclohexane", 0, {})
    taut = _taut("C1CCCCC1", "cyclohexane")
    assert parallel_check_nonarom_rings(taut, contnr) is taut


def test_parallel_check_nonarom_rings_discards_broken_ring() -> None:
    contnr = MolContainer("C1CCCCC1", "cyclohexane", 0, {})
    taut = _taut("c1ccccc1", "cyclohexane")
    assert parallel_check_nonarom_rings(taut, contnr) is None


def test_parallel_check_chiral_centers_keeps_matching_count() -> None:
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    taut = _taut("C[C@H](N)C(=O)O", "alanine")
    assert parallel_check_chiral_centers(taut, contnr) is taut


def test_parallel_check_chiral_centers_discards_changed_count() -> None:
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    taut = _taut("CCC(=O)O", "alanine")
    assert parallel_check_chiral_centers(taut, contnr) is None


def test_parallel_check_carbon_hydrogens_keeps_matching_count() -> None:
    contnr = MolContainer("CC(=O)C", "acetone", 0, {})
    taut = _taut("CC(=O)C", "acetone")
    assert parallel_check_carbon_hydrogens(taut, contnr) is taut


def test_parallel_check_carbon_hydrogens_discards_changed_count() -> None:
    contnr = MolContainer("CC(=O)C", "acetone", 0, {})
    taut = _taut("CC(O)=C", "acetone")
    assert parallel_check_carbon_hydrogens(taut, contnr) is None


def test_tauts_no_break_arom_rngs_filters_in_process() -> None:
    contnr = MolContainer("C1CCCCC1", "cyclohexane", 0, {})
    keep = _taut("C1CCCCC1", "cyclohexane")
    drop = _taut("c1ccccc1", "cyclohexane")
    result = tauts_no_break_arom_rngs([contnr], [keep, drop], 1, "serial", None)
    assert result == [keep]


def test_tauts_no_elim_chiral_filters_in_process() -> None:
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    keep = _taut("C[C@H](N)C(=O)O", "alanine")
    drop = _taut("CCC(=O)O", "alanine")
    result = tauts_no_elim_chiral([contnr], [keep, drop], 1, "serial", None)
    assert result == [keep]


def test_tauts_no_change_hs_to_cs_filters_in_process() -> None:
    contnr = MolContainer("CC(=O)C", "acetone", 0, {})
    keep = _taut("CC(=O)C", "acetone")
    drop = _taut("CC(O)=C", "acetone")
    result = tauts_no_change_hs_to_cs_unless_alpha_to_carbnyl(
        [contnr], [keep, drop], 1, "serial", None
    )
    assert result == [keep]