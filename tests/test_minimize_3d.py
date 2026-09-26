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
