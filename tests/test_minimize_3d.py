"""Unit tests for the final 3D minimization step."""

from gypsum_dl.steps.conf import Minimize3D


class _FakeConf:
    """A stand-in for MyConformer whose energy changes on minimize()."""

    def __init__(self, pre: float, post: float) -> None:
        self.energy = pre
        self._post = post
        self.minimized = False

    def minimize(self) -> None:
        if not self.minimized:
            self.energy = self._post
            self.minimized = True

    def conformer(self):
        # The energy carried on this token is what the (patched) MyConformer
        # reads back, so returning self is enough for the test.
        return self


class _FakeMol:
    def __init__(self, confs) -> None:
        self.conformers = confs
        self.genealogy = []

    def add_conformers(self, num, rmsd_cutoff=0.1, minimize=True) -> None:
        # Mirror the real add_conformers(minimize=False): it sorts by the
        # (pre-minimization) energy and does not minimize.
        self.conformers.sort(key=lambda c: c.energy)

    def smiles(self, noh=False):
        return "CCO"


def test_parallel_minit_returns_lowest_energy_minimized_conformer(monkeypatch) -> None:
    # Regression (M7): conformers are ranked by pre-minimization energy, then
    # the best few are minimized. Post-minimization ranking is not monotonic
    # with the pre-minimization ranking, so the code must re-sort the minimized
    # subset instead of blindly keeping conformers[0].
    #
    # Pre-min order: A(10) < B(20).  Post-min: B(1) < A(5). The lowest-energy
    # minimized conformer is B, but the buggy code returned A.
    conf_a = _FakeConf(pre=10, post=5)
    conf_b = _FakeConf(pre=20, post=1)
    mol = _FakeMol([conf_b, conf_a])  # deliberately unsorted

    captured = {}

    class _FakeMyConformer:
        def __init__(self, new_mol, conf, second_embed) -> None:
            self.energy = conf.energy
            captured["conf"] = conf

    monkeypatch.setattr(Minimize3D, "MyConformer", _FakeMyConformer)

    result = Minimize3D.parallel_minit(
        mol, max_variants_per_compound=2, thoroughness=1, second_embed=False
    )

    assert captured["conf"] is conf_b
    assert result.conformers[0].energy == 1
    assert "1 kcal/mol" in result.genealogy[-1]


class _AlertMol:
    """Minimal MyMol stand-in for the error-alert loop of minimize_3d."""

    def __init__(self, rdkit_mol) -> None:
        self.rdkit_mol = rdkit_mol
        self.contnr_idx = 0
        self.genealogy: list = []
        self.conformers: list = [object()]


class _AlertContnr:
    def __init__(self, mols) -> None:
        self.mols = mols
        # Non-zero so minimize_3d skips these mols (already minimized elsewhere)
        # and only its final error-alert loop runs — no RDKit work needed.
        self.num_nonaro_rngs = 1


def test_minimize_3d_flags_mols_that_failed_to_embed() -> None:
    # Regression (M2): a mol that failed 3D optimization has rdkit_mol == None,
    # never == "". The old `mol.rdkit_mol == ""` check never matched, so the
    # "Could not optimize 3D geometry" note and conformer clearing were dead
    # code. A None mol must now be flagged; a real mol must be left alone.
    failed = _AlertMol(None)
    ok = _AlertMol(object())
    contnr = _AlertContnr([failed, ok])

    Minimize3D.minimize_3d(
        [contnr],
        max_variants_per_compound=1,
        thoroughness=1,
        num_procs=1,
        second_embed=False,
        job_manager="serial",
        parallelizer_obj=None,
    )

    assert failed.genealogy[-1] == "(WARNING: Could not optimize 3D geometry)"
    assert failed.conformers == []
    assert ok.genealogy == []
    assert ok.conformers != []
