"""Unit tests for the final 3D minimization step."""

from typing import TypedDict

from gypsum_dl.steps.conf import Minimize3D, PrepareThreeD


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


_PrepareThreeDParams = TypedDict(
    "_PrepareThreeDParams",
    {
        "2d_output_only": bool,
        "max_variants_per_compound": int,
        "thoroughness": int,
        "num_processors": int,
        "job_manager": str,
        "Parallelizer": object | None,
        "second_embed": bool,
        "skip_alternate_ring_conformations": bool,
        "skip_optimize_geometry": bool,
    },
)


def _prepare_3d_params(skip_ring_confs: bool) -> _PrepareThreeDParams:
    """Build the parameter dict prepare_3d reads, for one ring-conf setting.

    Only the keys prepare_3d touches are included, so the test fails loudly if
    the step starts reading something new.

    Args:
        skip_ring_confs: Value for skip_alternate_ring_conformations.

    Returns:
        A parameter dict suitable for prepare_3d.
    """
    return {
        "2d_output_only": False,
        "max_variants_per_compound": 1,
        "thoroughness": 1,
        "num_processors": 1,
        "job_manager": "serial",
        "Parallelizer": None,
        "second_embed": False,
        "skip_alternate_ring_conformations": skip_ring_confs,
        "skip_optimize_geometry": False,
    }


def _stub_prepare_3d_steps(monkeypatch) -> dict[str, object]:
    """Replace the three steps prepare_3d calls with recording stubs.

    Lets the test check how prepare_3d wires its flags together without doing
    any RDKit work.

    Args:
        monkeypatch: The pytest monkeypatch fixture.

    Returns:
        A dict recording whether the ring-conformer step ran and what
        include_nonaro_rings value reached minimize_3d.
    """
    captured: dict[str, object] = {
        "ring_confs_ran": False,
        "include_nonaro_rings": None,
    }

    def fake_ring_confs(*args: object, **kwargs: object) -> None:
        captured["ring_confs_ran"] = True

    def fake_minimize_3d(*args: object, **kwargs: object) -> None:
        # The flag is accepted either way round, so the test does not break if
        # prepare_3d starts passing it positionally.
        positional = args[7] if len(args) > 7 else None
        captured["include_nonaro_rings"] = kwargs.get(
            "include_nonaro_rings", positional
        )

    monkeypatch.setattr(PrepareThreeD, "convert_2d_to_3d", lambda *a, **k: None)
    monkeypatch.setattr(
        PrepareThreeD, "generate_alternate_3d_nonaromatic_ring_confs", fake_ring_confs
    )
    monkeypatch.setattr(PrepareThreeD, "minimize_3d", fake_minimize_3d)
    return captured


def test_prepare_3d_minimizes_ring_mols_when_ring_conf_step_is_skipped(
    monkeypatch,
) -> None:
    # Regression: minimize_3d skips containers with non-aromatic rings, since
    # generating alternate ring conformations minimizes them as a side effect.
    # With --skip_alternate_ring_conformations, that side effect never happens,
    # so those molecules were written out with raw embedded coordinates and no
    # Energy property. prepare_3d must tell minimize_3d to handle them.
    captured = _stub_prepare_3d_steps(monkeypatch)

    PrepareThreeD.prepare_3d([], _prepare_3d_params(skip_ring_confs=True))

    assert captured["ring_confs_ran"] is False
    assert captured["include_nonaro_rings"] is True


def test_prepare_3d_leaves_ring_mols_to_the_ring_conf_step_by_default(
    monkeypatch,
) -> None:
    # The complement: when ring conformers are generated, those molecules are
    # already minimized, so minimize_3d must not redo them.
    captured = _stub_prepare_3d_steps(monkeypatch)

    PrepareThreeD.prepare_3d([], _prepare_3d_params(skip_ring_confs=False))

    assert captured["ring_confs_ran"] is True
    assert captured["include_nonaro_rings"] is False
