"""Regression tests for cis/trans double-bond enumeration."""

import random
import time

import pytest

from rdkit import Chem

from gypsum_dl.MyMol import MyMol
from gypsum_dl.steps.smiles.EnumerateDoubleBonds import (
    _could_carry_stereochemistry,
    parallel_get_double_bonded,
    sample_bond_dir_configs,
)


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


def test_preserves_the_input_smiles() -> None:
    # Regression (F2): the step set only contnr_idx, name, and genealogy on
    # each variant, so MyMol.__init__'s assumption that orig_smi is the
    # molecule's own SMILES stood. The PDB header's "Original SMILES string"
    # then repeated the variant instead of naming the library entry, and the
    # metal check in DurrantLabFilter reads orig_smi_deslt.
    mol = MyMol("CC=CC")
    mol.orig_smi = "CC=CC.[Na+]"
    mol.orig_smi_deslt = "CC=CC"
    mol.contnr_idx = 3
    mol.name = "butene"

    results = parallel_get_double_bonded(mol, 1, 1)

    assert len(results) == 2
    for result in results:
        assert result is not mol  # variants were actually built
        assert result.orig_smi == "CC=CC.[Na+]"
        assert result.orig_smi_deslt == "CC=CC"
        assert result.contnr_idx == 3
        assert result.name == "butene"


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


def test_output_does_not_carry_explicit_hydrogens() -> None:
    # Regression: enumeration works on Chem.AddHs(mol) and the candidate SMILES
    # were canonicalized from that molecule, so every variant came back spelled
    # with [H] atoms. MyMol parses with sanitize=False, which keeps those atoms
    # in the graph and in can_smi, the key used to tell variants apart, so the
    # same molecule spelled two ways consumed two max_variants_per_compound
    # slots. Tetrasubstituted alkene, so no hydrogen is needed to define the
    # stereochemistry and all of them should be gone.
    results = parallel_get_double_bonded(MyMol("CC(Cl)=C(Cl)C"), 1, 1)

    assert len(results) == 2
    for mol in results:
        assert "[H]" not in mol.smiles()
        assert all(atom.GetAtomicNum() != 1 for atom in mol.rdkit_mol.GetAtoms())


def test_carbonyls_do_not_crowd_out_a_stereogenic_bond() -> None:
    # Regression: the candidate list came from GetStereo() is STEREONONE, which
    # is also true of carbonyls, so four ketones competed with the one alkene
    # for a budget of a single bond. Four times out of five the shuffle kept a
    # carbonyl, no direction could be assigned, and the molecule came back with
    # its alkene geometry unspecified instead of enumerated.
    smi = "CC(=O)CC(=O)CC(=O)CC(=O)CC=CC"
    for seed in range(10):
        random.seed(seed)
        results = parallel_get_double_bonded(MyMol(smi), 1, 1)
        smis = [m.smiles(True) for m in results]
        assert len(results) == 2, f"seed {seed} produced {smis}"
        assert all("/" in s or "\\" in s for s in smis)


def test_log_counts_only_bonds_that_can_be_enumerated(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # The same miscount reached the log: the reported number of double bonds
    # "with unspecified stereochemistry" included the carbonyls, so a compound
    # with one alkene and four ketones claimed five.
    parallel_get_double_bonded(MyMol("CC(=O)CC(=O)CC(=O)CC(=O)CC=CC"), 1, 1)

    out = " ".join(capsys.readouterr().out.split())
    assert "has 1 double bond(s) with unspecified stereochemistry" in out
    assert "left unspecified" not in out


def test_terminal_alkene_does_not_crowd_out_a_stereogenic_bond() -> None:
    # Regression: the stereogenicity check counted bonds on the hydrogen-added
    # molecule, so a terminal =CH2 (two hydrogens) passed, competed with the
    # real alkene for a budget of one bond, and about half the seeds came back
    # with nothing enumerated and no warning.
    smi = "CC(=C)CC=CC"
    for seed in range(10):
        random.seed(seed)
        results = parallel_get_double_bonded(MyMol(smi), 1, 1)
        smis = [m.smiles(True) for m in results]
        assert len(results) == 2, f"seed {seed} produced {smis}"
        assert all("/" in s or "\\" in s for s in smis)


def test_terminal_alkene_is_not_logged_as_unspecified(
    capsys: pytest.CaptureFixture[str],
) -> None:
    mol = MyMol("CC=C")

    assert parallel_get_double_bonded(mol, 1, 1) == [mol]
    assert "unspecified stereochemistry" not in capsys.readouterr().out


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        ("CC=CC", True),
        ("CC=N", True),  # N-H imines are E/Z; a single hydrogen is enough.
        ("CC=C", False),
        ("CC(C)=C", False),
        ("CC=O", False),
    ],
)
def test_could_carry_stereochemistry(smiles: str, expected: bool) -> None:
    mol_with_hs = Chem.AddHs(Chem.MolFromSmiles(smiles))
    double_bond_idx = next(
        bond.GetIdx()
        for bond in mol_with_hs.GetBonds()
        if bond.GetBondType() == Chem.BondType.DOUBLE
    )

    assert _could_carry_stereochemistry(mol_with_hs, double_bond_idx) is expected


def test_sample_bond_dir_configs_enumerates_small_spaces() -> None:
    configs = sample_bond_dir_configs(4, 1024)

    assert len(set(configs)) == 16
    assert all(len(config) == 4 for config in configs)


def test_sample_bond_dir_configs_samples_large_spaces() -> None:
    # Regression: the direction space is 2**num_bonds, and num_bonds runs about
    # four times the number of double bonds being varied, so full enumeration
    # was exponential in a bond count that is only logarithmic in the variant
    # budget. Above the cap the space has to be sampled instead.
    configs = sample_bond_dir_configs(20, 500)

    assert len(set(configs)) == 500
    assert all(len(config) == 20 for config in configs)


def test_sample_bond_dir_configs_reproducible_when_seeded() -> None:
    random.seed(4321)
    first = sample_bond_dir_configs(20, 100)

    random.seed(4321)
    second = sample_bond_dir_configs(20, 100)

    assert first == second


def test_reports_when_only_some_double_bonds_are_varied(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Regression: num_bonds_to_keep truncates the unspecified double bonds to
    # a logarithmic count, at random, and the rest are never assigned a
    # direction (the 3D embedder picks their geometry arbitrarily later). The
    # log reported only the pre-truncation count, so a user reading "has 8
    # double bond(s) with unspecified stereochemistry" alongside a genealogy
    # entry of "(cis-trans isomerization)" had no way to learn that half of
    # them were left unspecified.
    smi = "CC=C" * 8 + "C"

    parallel_get_double_bonded(MyMol(smi), 5, 3)

    # Collapse whitespace: log() wraps at 80 columns.
    out = " ".join(capsys.readouterr().out.split())
    assert "only 4 of the 8 double bonds" in out
    assert "left unspecified" in out


def test_no_truncation_notice_when_all_bonds_are_varied(
    capsys: pytest.CaptureFixture[str],
) -> None:
    parallel_get_double_bonded(MyMol("CC=CC"), 1, 1)

    out = " ".join(capsys.readouterr().out.split())
    assert "has 1 double bond(s) with unspecified stereochemistry" in out
    assert "left unspecified" not in out


def test_discards_variants_that_cannot_be_canonicalized(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: the guard read new_mol.can_smi, which is still "" for a
    # MyMol built from a SMILES string (only a later smiles() call can set it
    # to None), so the condition was always true and the guard was dead. A
    # variant whose canonicalization fails then reached the results carrying
    # can_smi = None, which leaked into the genealogy strings and into the
    # set-based deduplication downstream.
    original_smiles = MyMol.smiles

    def failing_smiles(self: MyMol, noh: bool = False) -> str | None:
        """Fail canonicalization for the enumerated variants only.

        The variants are the molecules built from stereochemistry-assigned
        SMILES, so they are the ones whose SMILES carry a slash or backslash.
        The input molecule has to keep working, because it is canonicalized
        for the log message before the guard is ever reached.

        Args:
            self: The molecule being canonicalized.
            noh: Whether to omit hydrogens, as in `MyMol.smiles`.

        Returns:
            None for an enumerated variant, otherwise the canonical SMILES.
        """
        if not noh and ("/" in self.orig_smi or "\\" in self.orig_smi):
            return None
        return original_smiles(self, noh)

    assert parallel_get_double_bonded(MyMol("CC=CC"), 1, 1)

    monkeypatch.setattr(MyMol, "smiles", failing_smiles)

    assert parallel_get_double_bonded(MyMol("CC=CC"), 1, 1) == []


def test_many_unspecified_double_bonds_stays_bounded() -> None:
    # Regression: at these settings num_bonds_to_keep is 5, each retained double
    # bond contributes up to four single bonds, and the full product of 2**20
    # direction assignments was materialized as a list and then stereo-assigned
    # one at a time. The call took minutes and hundreds of megabytes.
    smi = "CC=C" * 8 + "C"

    started = time.monotonic()
    results = parallel_get_double_bonded(MyMol(smi), 16, 2)
    elapsed = time.monotonic() - started

    assert results
    assert elapsed < 30
