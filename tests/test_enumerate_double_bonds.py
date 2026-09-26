"""Regression tests for cis/trans double-bond enumeration."""

from gypsum_dl.MyMol import MyMol
from gypsum_dl.steps.smiles.EnumerateDoubleBonds import parallel_get_double_bonded


def test_enumerates_when_thoroughness_times_variants_is_one() -> None:
    # Regression: with thoroughness * max_variants_per_compound == 1, log2(1) == 0
    # used to zero out num_bonds_to_keep, silently disabling enumeration.
    results = parallel_get_double_bonded(MyMol("CC=CC"), 1, 1)
    smis = [m.smiles(True) for m in results]
    assert len(results) == 2
    assert all("/" in s or "\\" in s for s in smis)
