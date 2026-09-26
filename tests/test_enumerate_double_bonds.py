"""Regression tests for cis/trans double-bond enumeration."""

import random

from gypsum_dl.MyMol import MyMol
from gypsum_dl.steps.smiles.EnumerateDoubleBonds import parallel_get_double_bonded


def test_enumerates_when_thoroughness_times_variants_is_one() -> None:
    # Regression: with thoroughness * max_variants_per_compound == 1, log2(1) == 0
    # used to zero out num_bonds_to_keep, silently disabling enumeration.
    results = parallel_get_double_bonded(MyMol("CC=CC"), 1, 1)
    smis = [m.smiles(True) for m in results]
    assert len(results) == 2
    assert all("/" in s or "\\" in s for s in smis)


def test_does_not_mutate_input_molecule() -> None:
    # Regression (M9): parallel_get_double_bonded used to do
    # `mol.rdkit_mol = Chem.AddHs(mol.rdkit_mol)`, mutating a molecule owned by
    # the container. In serial mode the extra explicit hydrogens leaked into
    # later steps; under multiprocessing they did not, so serial and
    # multiprocessing diverged. The molecule must be left untouched.
    mol = MyMol("CC=CC")
    n_atoms_before = mol.rdkit_mol.GetNumAtoms()

    parallel_get_double_bonded(mol, 1, 1)

    assert mol.rdkit_mol.GetNumAtoms() == n_atoms_before


def test_reproducible_when_seeded() -> None:
    # Regression (M10): selection of which double bonds to enumerate uses
    # random.shuffle, and the candidate SMILES were pulled from a set whose
    # iteration order varies with PYTHONHASHSEED. Given the same seed, the
    # result must be identical run to run. Uses more unspecified double bonds
    # (5) than num_bonds_to_keep (4 for these args) so the shuffle actually
    # chooses a subset.
    smi = "CC=CCC=CCC=CCC=CCC=CC"

    random.seed(1234)
    res1 = sorted(m.smiles(True) for m in parallel_get_double_bonded(MyMol(smi), 5, 3))

    random.seed(1234)
    res2 = sorted(m.smiles(True) for m in parallel_get_double_bonded(MyMol(smi), 5, 3))

    assert res1  # enumeration actually produced something
    assert res1 == res2
