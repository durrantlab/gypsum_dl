"""Regression tests for chiral enumeration."""

import pytest

from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.MyMol import MyMol
from gypsum_dl.steps.smiles import EnumerateChiralMols
from gypsum_dl.steps.smiles.EnumerateChiralMols import (
    enumerate_chiral_molecules,
    parallel_get_chiral,
)


def test_parallel_get_chiral_receives_args_in_declared_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression (bug 9): the call site built params as
    # (mol, thoroughness, max_variants_per_compound) while parallel_get_chiral
    # is declared (mol, max_variants_per_compound, thoroughness). The product of
    # the two hid it; distinct values expose the swap.
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    contnr.add_smiles("C[C@H](N)C(=O)O")
    mol = contnr.mols[0]

    captured = []

    def spy(m, max_variants_per_compound, thoroughness):
        captured.append((m, max_variants_per_compound, thoroughness))
        return [m]

    monkeypatch.setattr(EnumerateChiralMols, "parallel_get_chiral", spy)

    enumerate_chiral_molecules([contnr], 3, 5, 1, "serial", None)

    assert captured == [(mol, 3, 5)]


MANY_UNSPECIFIED_CENTERS = "CC(N)C(O)C(F)C(Cl)C(Br)C"


def _mol_with_many_unspecified_centers() -> tuple[MyMol, int]:
    """Build a molecule whose unspecified chiral centers exceed any small cap.

    Returns:
        The molecule and the number of its unspecified chiral centers, so the
        tests can assert against the count RDKit actually perceives rather than
        a hard-coded one.
    """
    contnr = MolContainer(MANY_UNSPECIFIED_CENTERS, "polychiral", 0, {})
    contnr.add_smiles(MANY_UNSPECIFIED_CENTERS)
    mol = contnr.mols[0]
    num = len([p for p in mol.chiral_cntrs_w_unasignd() if p[1] == "?"])
    return mol, num


def test_parallel_get_chiral_completes_truncated_assignments(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: the combinatorial cap broke out of the assignment loop early,
    # so each option covered only the leading centers. zip() then dropped the
    # rest, emitting variants with unspecified chirality.
    mol, num = _mol_with_many_unspecified_centers()
    assert num > 1

    captured: list[list[list[str]]] = []

    def spy(lst: list[list[str]], keep: int, msg_if_cut: str = "") -> list[list[str]]:
        captured.append([list(option) for option in lst])
        return lst[:keep]

    monkeypatch.setattr(EnumerateChiralMols.utils, "random_sample", spy)

    parallel_get_chiral(mol, 1, 1)

    assert len(captured) == 1
    assert captured[0]
    assert all(len(option) == num for option in captured[0])


def test_parallel_get_chiral_preserves_the_input_smiles() -> None:
    # Regression (F2): the step set only contnr_idx, name, and genealogy on
    # each variant, so MyMol.__init__'s assumption that orig_smi is the
    # molecule's own SMILES stood. The PDB header's "Original SMILES string"
    # then repeated the variant instead of naming the library entry, and the
    # metal check in DurrantLabFilter reads orig_smi_deslt.
    mol = MyMol("CC(N)C(=O)O")
    mol.orig_smi = "CC(N)C(=O)O.[Na+]"
    mol.orig_smi_deslt = "CC(N)C(=O)O"
    mol.contnr_idx = 3
    mol.name = "alanine"

    results = parallel_get_chiral(mol, 2, 1)

    assert results
    for result in results:
        assert result is not mol  # variants were actually built
        assert result.orig_smi == "CC(N)C(=O)O.[Na+]"
        assert result.orig_smi_deslt == "CC(N)C(=O)O"
        assert result.contnr_idx == 3
        assert result.name == "alanine"


def test_parallel_get_chiral_emits_fully_specified_variants() -> None:
    mol, num = _mol_with_many_unspecified_centers()
    assert num > 1

    results = parallel_get_chiral(mol, 1, 1)

    assert results
    for result in results:
        unassigned = [p for p in result.chiral_cntrs_w_unasignd() if p[1] == "?"]
        assert unassigned == []
