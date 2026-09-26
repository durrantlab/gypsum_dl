"""Regression tests for chiral enumeration."""

import pytest

from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.steps.smiles import EnumerateChiralMols
from gypsum_dl.steps.smiles.EnumerateChiralMols import enumerate_chiral_molecules


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
