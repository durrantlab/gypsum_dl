"""Tests for the tautomer rejection filters."""

import pytest

from gypsum_dl import MyMol
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.steps.smiles import MakeTautomers
from gypsum_dl.steps.smiles.MakeTautomers import (
    parallel_check_chiral_centers,
    parallel_check_nonarom_rings,
    parallel_make_taut,
    tauts_no_break_arom_rngs,
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


def test_parallel_check_nonarom_rings_discards_changed_aromaticity() -> None:
    # The criterion is symmetric: aromatizing a ring that was nonaromatic is
    # rejected just as dearomatizing an aromatic one is.
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


def test_tauts_no_break_arom_rngs_drops_orphan_taut() -> None:
    # Regression (bug 10): a taut matching no container silently reused the last
    # (stale) `container`, so it could be compared against the wrong molecule
    # and kept. It must be dropped instead.
    contnr = MolContainer("c1ccccc1", "benzene", 0, {})
    orphan = _taut("c1ccccc1", "benzene")
    orphan.contnr_idx = 99
    result = tauts_no_break_arom_rngs([contnr], [orphan], 1, "serial", None)
    assert result == []


def test_tauts_no_elim_chiral_drops_orphan_taut() -> None:
    # Regression (bug 10): same stale-container reuse in the chiral filter.
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    orphan = _taut("C[C@H](N)C(=O)O", "alanine")
    orphan.contnr_idx = 99
    result = tauts_no_elim_chiral([contnr], [orphan], 1, "serial", None)
    assert result == []


def test_enol_tautomers_survive_the_filters() -> None:
    # A retired filter compared the total count of hydrogens bound to carbon
    # and rejected any change, which discards every keto-enol pair (acetone
    # carries six such hydrogens, its enol five). Keto-enol tautomerism is
    # documented behavior, so the surviving filters must let the enol through.
    contnr = MolContainer("CC(=O)C", "acetone", 0, {})
    keto = _taut("CC(=O)C", "acetone")
    enol = _taut("CC(O)=C", "acetone")

    kept = tauts_no_break_arom_rngs([contnr], [keto, enol], 1, "serial", None)
    kept = tauts_no_elim_chiral([contnr], kept, 1, "serial", None)

    assert kept == [keto, enol]
