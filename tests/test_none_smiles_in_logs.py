"""Regression tests for log messages that concatenate MyMol.smiles(True).

smiles(True) returns None when a variant cannot be deprotonated or
canonicalized. Messages built with bare string concatenation raised TypeError
on that None, and inside a parallelized step the exception discarded the
whole container rather than logging the one variant.
"""

from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.MyMol import MyMol
from gypsum_dl.steps.smiles.DurrantLabFilter import parallel_durrant_lab_filter
from gypsum_dl.steps.smiles.EnumerateChiralMols import parallel_get_chiral
from gypsum_dl.steps.smiles.EnumerateDoubleBonds import parallel_get_double_bonded
from gypsum_dl.steps.smiles.MakeTautomers import (
    chirality_facts,
    parallel_check_chiral_centers,
    parallel_check_nonarom_rings,
    parallel_make_taut,
    ring_facts,
)


def _variant_without_noh_smiles(smiles: str, name: str) -> MyMol:
    """Build a variant whose hydrogen-free SMILES has already failed.

    Caching None is what smiles(True) itself does when deprotonation fails,
    so this reproduces that state without depending on which inputs happen to
    break RemoveHs in the installed RDKit.

    Args:
        smiles: SMILES string for the variant.
        name: Ligand name used in the log messages.

    Returns:
        A MyMol whose smiles(True) returns None.
    """
    mol = MyMol(smiles)
    mol.name = name
    mol.contnr_idx = 0
    mol.can_smi_noh = None
    assert mol.rdkit_mol is not None
    assert mol.smiles(True) is None
    return mol


def test_durrant_filter_keeps_good_variant_beside_one_without_smiles() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    good = MyMol("CCO")
    good.name = "ethanol"
    good.contnr_idx = 0
    # Boron trips the [#5] pattern, so this variant reaches the rejection log.
    bad = _variant_without_noh_smiles("CB(O)O", "ethanol")
    contnr.mols = [good, bad]

    result = parallel_durrant_lab_filter(contnr)

    assert result is contnr
    assert len(result.mols) == 1
    assert result.mols[0] is good


def test_chiral_enumeration_logs_variant_without_noh_smiles() -> None:
    mol = _variant_without_noh_smiles("CC(N)C(=O)O", "alanine")
    assert parallel_get_chiral(mol, 2, 1)


def test_double_bond_enumeration_logs_variant_without_noh_smiles() -> None:
    mol = _variant_without_noh_smiles("CC=CC", "butene")
    assert len(parallel_get_double_bonded(mol, 1, 1)) == 2


def test_make_taut_logs_variant_without_noh_smiles() -> None:
    contnr = MolContainer("Oc1ccccn1", "hydroxypyridine", 0, {})
    mol = _variant_without_noh_smiles("Oc1ccccn1", "hydroxypyridine")

    results = parallel_make_taut(mol, contnr.contnr_props(), 10)

    # More than one form is what routes through the "has tautomers" message.
    assert results is not None
    assert len(results) > 1


def test_nonarom_ring_rejection_logs_tautomer_without_noh_smiles() -> None:
    contnr = MolContainer("C1CCCCC1", "cyclohexane", 0, {})
    taut = _variant_without_noh_smiles("c1ccccc1", "cyclohexane")
    assert parallel_check_nonarom_rings(taut, ring_facts(contnr)) is None


def test_chiral_center_rejection_logs_tautomer_without_noh_smiles() -> None:
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    taut = _variant_without_noh_smiles("CCC(=O)O", "alanine")
    assert parallel_check_chiral_centers(taut, chirality_facts(contnr)) is None
