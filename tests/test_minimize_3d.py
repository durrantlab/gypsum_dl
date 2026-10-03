"""Unit tests for the final 3D minimization step."""

from typing import TypedDict

import pytest

from gypsum_dl import MyMol
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

    def add_conformers(
        self, num, rmsd_cutoff=0.1, minimize=True, second_embed=False
    ) -> None:
        # Mirror the real add_conformers(minimize=False): it sorts by the
        # (pre-minimization) energy and does not minimize.
        self.conformers.sort(key=lambda c: c.energy)

    def smiles(self, noh=False):
        return "CCO"


def test_parallel_minit_returns_lowest_energy_minimized_conformer(monkeypatch) -> None:
    # Regression (M7): conformers arrive sorted by pre-minimization energy.
    # Post-minimization ranking is not monotonic with the pre-minimization
    # ranking, so the code must re-sort the minimized conformers instead of
    # blindly keeping conformers[0].
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

    # parallel_minit works on a copy of mol, so B is recognized by the energy
    # only it reaches after minimization rather than by identity.
    assert captured["conf"].energy == 1
    assert result.conformers[0].energy == 1
    assert "1 kcal/mol" in result.genealogy[-1]


def test_parallel_minit_minimizes_conformers_beyond_the_variant_cap(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Thoroughness must buy minimized candidates, not just discarded embeds.

    Regression: only the top max_variants_per_compound conformers by raw
    (unminimized) energy were minimized. With a cap of one, the extra
    conformers thoroughness generated were thrown away unminimized, so a
    conformer that embeds strained but relaxes lowest could never win. Here C
    has the worst raw energy and the best minimized energy.
    """
    conf_a = _FakeConf(pre=10, post=5)
    conf_b = _FakeConf(pre=20, post=3)
    conf_c = _FakeConf(pre=30, post=1)
    mol = _FakeMol([conf_c, conf_a, conf_b])

    class _FakeMyConformer:
        def __init__(
            self, new_mol: object, conf: _FakeConf, second_embed: bool
        ) -> None:
            self.energy = conf.energy

    monkeypatch.setattr(Minimize3D, "MyConformer", _FakeMyConformer)

    result = Minimize3D.parallel_minit(
        mol, max_variants_per_compound=1, thoroughness=3, second_embed=False
    )

    assert result.conformers[0].energy == 1
    assert "1 kcal/mol" in result.genealogy[-1]


def test_parallel_minit_genealogy_omits_the_failure_sentinel(monkeypatch) -> None:
    # Regression: a conformer whose force field could not be set up carries an
    # infinite energy, so the genealogy line read "optimized conformer: inf
    # kcal/mol", which looks like a measurement.
    conf = _FakeConf(pre=float("inf"), post=float("inf"))
    mol = _FakeMol([conf])

    class _FakeMyConformer:
        def __init__(self, new_mol, conf, second_embed) -> None:
            self.energy = conf.energy

    monkeypatch.setattr(Minimize3D, "MyConformer", _FakeMyConformer)

    result = Minimize3D.parallel_minit(
        mol, max_variants_per_compound=1, thoroughness=1, second_embed=False
    )

    assert "inf" not in result.genealogy[-1]
    assert "force field failed" in result.genealogy[-1]


class _EnergyMol:
    """Minimal MyMol stand-in carrying one conformer energy and a props dict."""

    def __init__(self, energy: float) -> None:
        self.conformers: list = [_FakeConf(pre=energy, post=energy)]
        self.mol_props: dict = {}
        self.genealogy: list = []
        self.contnr_idx = 0
        self.rdkit_mol = object()


class _EnergyContnr:
    def __init__(self, mols: list) -> None:
        self.mols = mols
        self.contnr_idx = 0
        # Zero so minimize_3d handles these molecules rather than leaving them
        # to the ring-conformer step.
        self.num_nonaro_rngs = 0

    def add_mol(self, mol) -> None:
        self.mols.append(mol)


def test_minimize_3d_leaves_out_the_energy_property_when_scoring_failed(
    monkeypatch,
) -> None:
    # Regression: the infinite sentinel was written straight into the SDF
    # Energy field. A None value is skipped by set_rdkit_mol_prop, so the field
    # is absent instead of holding an unparseable number.
    failed = _EnergyMol(float("inf"))
    scored = _EnergyMol(-3.5)
    contnr = _EnergyContnr([failed, scored])

    monkeypatch.setattr(Minimize3D, "parallel_minit", lambda mol, *a, **k: mol)

    Minimize3D.minimize_3d(
        [contnr],
        max_variants_per_compound=1,
        thoroughness=1,
        num_procs=1,
        second_embed=False,
        job_manager="serial",
        parallelizer_obj=None,
    )

    assert failed.mol_props["Energy"] is None
    assert scored.mol_props["Energy"] == -3.5


class _RecordingMol(_FakeMol):
    """A _FakeMol that remembers how many conformers it was asked for.

    The zero-cap regression is about the count that reaches add_conformers (and
    the slice taken afterwards), so the test needs to see that number rather
    than only the returned molecule.
    """

    def __init__(self, confs: list) -> None:
        super().__init__(confs)
        self.requested = -1
        self.second_embed: bool | None = None

    def add_conformers(
        self, num, rmsd_cutoff=0.1, minimize=True, second_embed=False
    ) -> None:
        self.requested = num
        self.second_embed = second_embed
        super().add_conformers(num, rmsd_cutoff, minimize, second_embed)


def test_parallel_minit_keeps_one_conformer_with_a_zero_variant_cap(
    monkeypatch,
) -> None:
    # Regression: max_variants_per_compound is allowed to be 0, which the
    # SMILES enumeration steps read as "do not enumerate variants." Here it
    # requested zero conformers and then indexed conformers[:0][0], so the
    # IndexError was swallowed by the worker wrapper and the molecule was
    # reported as one that simply produced nothing.
    conf = _FakeConf(pre=10, post=5)
    mol = _RecordingMol([conf])

    class _FakeMyConformer:
        def __init__(self, new_mol, conf, second_embed) -> None:
            self.energy = conf.energy

    monkeypatch.setattr(Minimize3D, "MyConformer", _FakeMyConformer)

    result = Minimize3D.parallel_minit(
        mol, max_variants_per_compound=0, thoroughness=1, second_embed=False
    )

    # The request is recorded on the copy parallel_minit works on, which the
    # returned molecule is deep-copied from.
    assert result is not None
    assert result.requested >= 1
    assert len(result.conformers) == 1


def test_parallel_minit_forwards_second_embed(monkeypatch) -> None:
    # Regression: parallel_minit received second_embed but never handed it to
    # add_conformers, so the fallback embedder could not run from this step.
    mol = _RecordingMol([_FakeConf(pre=10, post=5)])

    class _FakeMyConformer:
        def __init__(self, new_mol, conf, second_embed) -> None:
            self.energy = conf.energy

    monkeypatch.setattr(Minimize3D, "MyConformer", _FakeMyConformer)

    result = Minimize3D.parallel_minit(
        mol, max_variants_per_compound=1, thoroughness=1, second_embed=True
    )

    assert result is not None
    assert result.second_embed is True


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
        # Every real container carries an index, and minimize_3d now maps the
        # containers by it rather than trusting list position.
        self.contnr_idx = 0
        # Non-zero so minimize_3d skips these mols (already minimized elsewhere)
        # and only its final error-alert loop runs: no RDKit work needed.
        self.num_nonaro_rngs = 1


def test_parallel_minit_leaves_its_input_untouched(monkeypatch) -> None:
    # Regression: serial and in-process runs hand parallel_minit the
    # container's own molecule. It sorted that molecule's conformers and
    # minimized them in place, and when every variant of a container failed,
    # minimize_3d kept those mutated originals, so the geometry written out
    # depended on the job manager.
    conf_a = _FakeConf(pre=10, post=5)
    conf_b = _FakeConf(pre=20, post=1)
    mol = _FakeMol([conf_b, conf_a])

    class _FakeMyConformer:
        def __init__(self, new_mol, conf, second_embed) -> None:
            self.energy = conf.energy

    monkeypatch.setattr(Minimize3D, "MyConformer", _FakeMyConformer)

    Minimize3D.parallel_minit(
        mol, max_variants_per_compound=2, thoroughness=1, second_embed=False
    )

    assert mol.conformers == [conf_b, conf_a]
    assert not conf_a.minimized
    assert not conf_b.minimized
    assert mol.genealogy == []


def test_minimize_3d_minimizes_containers_the_ring_step_could_not_process(
    monkeypatch,
) -> None:
    # Regression: minimize_3d skips ring-bearing containers because the
    # ring-conformer step minimizes them, but a container that step failed on
    # was never minimized anywhere and shipped its raw embedded geometry.
    contnr = _EnergyContnr([_EnergyMol(-3.5)])
    contnr.num_nonaro_rngs = 1
    minimized: list = []

    def fake_parallel_minit(mol, *args: object) -> object:
        minimized.append(mol)
        return mol

    monkeypatch.setattr(Minimize3D, "parallel_minit", fake_parallel_minit)

    Minimize3D.minimize_3d(
        [contnr],
        max_variants_per_compound=1,
        thoroughness=1,
        num_procs=1,
        second_embed=False,
        job_manager="serial",
        parallelizer_obj=None,
        ring_conf_failed_contnr_idxs=frozenset({0}),
    )

    assert minimized == contnr.mols
    assert contnr.mols[0].mol_props["Energy"] == -3.5


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

    # Regression: the flagged molecule also has to leave the container. It is
    # not writable (load_conformers_into_rdkit_mol returns early on a None
    # rdkit_mol), but the steps that follow minimization call accessors on it
    # from the main process, where an exception ends the run after all the
    # expensive work is done.
    assert contnr.mols == [ok]


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


def _prepare_3d_params(
    skip_ring_confs: bool, skip_optimize: bool = False
) -> _PrepareThreeDParams:
    """Build the parameter dict prepare_3d reads, for one ring-conf setting.

    Only the keys prepare_3d touches are included, so the test fails loudly if
    the step starts reading something new.

    Args:
        skip_ring_confs: Value for skip_alternate_ring_conformations.
        skip_optimize: Value for skip_optimize_geometry.

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
        "skip_optimize_geometry": skip_optimize,
    }


_STUB_RING_CONF_FAILURES: frozenset[int] = frozenset({7})


def _stub_prepare_3d_steps(monkeypatch) -> dict[str, object]:
    """Replace the three steps prepare_3d calls with recording stubs.

    Lets the test check how prepare_3d wires its flags together without doing
    any RDKit work.

    Args:
        monkeypatch: The pytest monkeypatch fixture.

    Returns:
        A dict recording whether the ring-conformer step ran, the minimize
        value it received, and what include_nonaro_rings and
        ring_conf_failed_contnr_idxs values reached minimize_3d. The stubbed
        ring-conformer step reports _STUB_RING_CONF_FAILURES as its failures.
    """
    captured: dict[str, object] = {
        "ring_confs_ran": False,
        "ring_confs_minimize": None,
        "include_nonaro_rings": None,
        "ring_conf_failed_contnr_idxs": None,
    }

    def fake_ring_confs(*args: object, **kwargs: object) -> frozenset[int]:
        captured["ring_confs_ran"] = True
        # Accepted either way round, so the test does not break if prepare_3d
        # starts passing the flag positionally.
        positional = args[7] if len(args) > 7 else None
        captured["ring_confs_minimize"] = kwargs.get("minimize", positional)
        return _STUB_RING_CONF_FAILURES

    def fake_minimize_3d(*args: object, **kwargs: object) -> None:
        # The flag is accepted either way round, so the test does not break if
        # prepare_3d starts passing it positionally.
        positional = args[7] if len(args) > 7 else None
        captured["include_nonaro_rings"] = kwargs.get(
            "include_nonaro_rings", positional
        )
        captured["ring_conf_failed_contnr_idxs"] = kwargs.get(
            "ring_conf_failed_contnr_idxs", args[8] if len(args) > 8 else None
        )

    monkeypatch.setattr(PrepareThreeD, "convert_2d_to_3d", lambda *a, **k: None)
    monkeypatch.setattr(
        PrepareThreeD, "generate_alternate_3d_nonaromatic_ring_confs", fake_ring_confs
    )
    monkeypatch.setattr(PrepareThreeD, "minimize_3d", fake_minimize_3d)
    return captured


def test_prepare_3d_does_not_minimize_ring_confs_when_optimization_is_skipped(
    monkeypatch,
) -> None:
    # Regression: the ring-conformer step calls add_conformers(..., minimize)
    # with a hardcoded True, so --skip_optimize_geometry skipped minimize_3d
    # but still returned UFF-minimized geometries for every molecule with a
    # non-aromatic ring, which is most drug-like input. prepare_3d has to pass
    # the flag along.
    captured = _stub_prepare_3d_steps(monkeypatch)

    PrepareThreeD.prepare_3d(
        [], _prepare_3d_params(skip_ring_confs=False, skip_optimize=True)
    )

    assert captured["ring_confs_ran"] is True
    assert captured["ring_confs_minimize"] is False
    # The separate minimization step is skipped outright, as before.
    assert captured["include_nonaro_rings"] is None


def test_prepare_3d_minimizes_ring_confs_by_default(monkeypatch) -> None:
    # The complement: without the flag, the ring conformers are still
    # minimized, so existing output does not change.
    captured = _stub_prepare_3d_steps(monkeypatch)

    PrepareThreeD.prepare_3d([], _prepare_3d_params(skip_ring_confs=False))

    assert captured["ring_confs_minimize"] is True


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


def test_prepare_3d_hands_ring_conf_failures_to_minimize_3d(monkeypatch) -> None:
    # Regression: containers the ring-conformer step could not process were
    # skipped by minimize_3d as well, so they were never minimized. prepare_3d
    # has to pass that step's failures along.
    captured = _stub_prepare_3d_steps(monkeypatch)

    PrepareThreeD.prepare_3d([], _prepare_3d_params(skip_ring_confs=False))

    assert captured["ring_conf_failed_contnr_idxs"] == _STUB_RING_CONF_FAILURES


def test_parallel_minit_survives_a_mol_with_no_rdkit_mol() -> None:
    # Regression (F6): the real MyConformer left .energy unset on its failure
    # paths, and parallel_minit reads it unconditionally. A molecule that
    # reached this step with conformers but a None rdkit_mol therefore raised
    # AttributeError inside the worker, which reported it as a molecule that
    # simply produced nothing. It must come back with the failed conformer and
    # a genealogy line saying so; minimize_3d's own alert loop then drops it.
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    assert mol.conformers
    mol.rdkit_mol = None

    result = Minimize3D.parallel_minit(
        mol, max_variants_per_compound=1, thoroughness=1, second_embed=False
    )

    assert result is not None
    assert result.conformers[0].mol is False
    assert result.conformers[0].energy == float("inf")
    assert "energy unavailable" in result.genealogy[-1]
