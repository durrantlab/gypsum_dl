"""Tests for the tautomer rejection filters."""

import pytest

from gypsum_dl import MyMol
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.steps.smiles import MakeTautomers
from gypsum_dl.steps.smiles.MakeTautomers import (
    chirality_facts,
    parallel_check_chiral_centers,
    parallel_check_nonarom_rings,
    parallel_make_taut,
    ring_facts,
    tauts_no_break_arom_rngs,
    tauts_no_elim_chiral,
)


def _taut(smiles: str, name: str) -> MyMol.MyMol:
    """Build a candidate tautomer tagged for container zero.

    Args:
        smiles: SMILES string for the tautomer.
        name: Ligand name used in the rejection log messages.

    Returns:
        The tagged MyMol instance.
    """
    mol = MyMol.MyMol(smiles)
    mol.name = name
    mol.contnr_idx = 0
    return mol


def test_parallel_make_taut_returns_none_when_rdkit_mol_unsanitizable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: parallel_make_taut called Chem.RemoveHs before its None
    # check, so an unsanitizable molecule (rdkit_mol is None) raised inside
    # RemoveHs instead of being dropped gracefully.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.mols.append(_taut("CCO", "ethanol"))

    class _NoneMol:
        rdkit_mol = None

    monkeypatch.setattr(MakeTautomers.MyMol, "MyMol", lambda *a, **k: _NoneMol())
    assert parallel_make_taut(contnr.mols[0], contnr.contnr_props(), 1) is None


def test_parallel_check_nonarom_rings_keeps_matching_tautomer() -> None:
    contnr = MolContainer("C1CCCCC1", "cyclohexane", 0, {})
    taut = _taut("C1CCCCC1", "cyclohexane")
    assert parallel_check_nonarom_rings(taut, ring_facts(contnr)) is taut


def test_parallel_check_nonarom_rings_discards_changed_aromaticity() -> None:
    # The criterion is symmetric: aromatizing a ring that was nonaromatic is
    # rejected just as dearomatizing an aromatic one is.
    contnr = MolContainer("C1CCCCC1", "cyclohexane", 0, {})
    taut = _taut("c1ccccc1", "cyclohexane")
    assert parallel_check_nonarom_rings(taut, ring_facts(contnr)) is None


def test_parallel_check_chiral_centers_keeps_matching_count() -> None:
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    taut = _taut("C[C@H](N)C(=O)O", "alanine")
    assert parallel_check_chiral_centers(taut, chirality_facts(contnr)) is taut


def test_parallel_check_chiral_centers_discards_changed_count() -> None:
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    taut = _taut("CCC(=O)O", "alanine")
    assert parallel_check_chiral_centers(taut, chirality_facts(contnr)) is None


def test_tauts_no_break_arom_rngs_filters_in_process() -> None:
    contnr = MolContainer("C1CCCCC1", "cyclohexane", 0, {})
    keep = _taut("C1CCCCC1", "cyclohexane")
    drop = _taut("c1ccccc1", "cyclohexane")
    result = tauts_no_break_arom_rngs([contnr], [keep, drop], 1, "serial", None)
    assert result == [keep]


def test_tauts_no_elim_chiral_filters_in_process() -> None:
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    keep = _taut("C[C@H](N)C(=O)O", "alanine")
    drop = _taut("CCC(=O)O", "alanine")
    result = tauts_no_elim_chiral([contnr], [keep, drop], 1, "serial", None)
    assert result == [keep]


def test_tauts_no_break_arom_rngs_drops_orphan_taut() -> None:
    # Regression (bug 10): a taut matching no container silently reused the last
    # (stale) `container`, so it could be compared against the wrong molecule
    # and kept. It must be dropped instead.
    contnr = MolContainer("c1ccccc1", "benzene", 0, {})
    orphan = _taut("c1ccccc1", "benzene")
    orphan.contnr_idx = 99
    result = tauts_no_break_arom_rngs([contnr], [orphan], 1, "serial", None)
    assert result == []


def test_tauts_no_elim_chiral_drops_orphan_taut() -> None:
    # Regression (bug 10): same stale-container reuse in the chiral filter.
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    orphan = _taut("C[C@H](N)C(=O)O", "alanine")
    orphan.contnr_idx = 99
    result = tauts_no_elim_chiral([contnr], [orphan], 1, "serial", None)
    assert result == []


def test_enol_tautomers_survive_the_filters() -> None:
    # A retired filter compared the total count of hydrogens bound to carbon
    # and rejected any change, which discards every keto-enol pair (acetone
    # carries six such hydrogens, its enol five). Keto-enol tautomerism is
    # documented behavior, so the surviving filters must let the enol through.
    contnr = MolContainer("CC(=O)C", "acetone", 0, {})
    keto = _taut("CC(=O)C", "acetone")
    enol = _taut("CC(O)=C", "acetone")

    kept = tauts_no_break_arom_rngs([contnr], [keto, enol], 1, "serial", None)
    kept = tauts_no_elim_chiral([contnr], kept, 1, "serial", None)

    assert kept == [keto, enol]


def _logged_message(captured: str) -> str:
    """Collapse a captured log record into one space-separated line.

    utils.log wraps at 80 columns, so an assertion on message wording has to
    ignore where the line breaks landed.

    Args:
        captured: Everything the log sink wrote during the call under test.

    Returns:
        The same text with every run of whitespace reduced to one space.
    """
    return " ".join(captured.split())


def test_nonarom_ring_rejection_names_non_aromatic_rings(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Regression (F7): the filter compares counts of non-aromatic rings, but
    # the rejection message reported a change in the number of aromatic rings,
    # sending anyone chasing a missing tautomer after the wrong quantity.
    contnr = MolContainer("C1CCCCC1", "cyclohexane", 0, {})
    taut = _taut("c1ccccc1", "cyclohexane")

    assert parallel_check_nonarom_rings(taut, ring_facts(contnr)) is None

    message = _logged_message(capsys.readouterr().out)
    assert "changed the number of non-aromatic rings" in message


def test_chiral_rejection_does_not_call_the_count_specified(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Regression (F7): the message called the number it printed the total
    # number of specified centers, which was neither of the things being
    # compared. Also a regression on the arithmetic behind that number: the
    # filter summed chiral_cntrs_only_asignd() and chiral_cntrs_w_unasignd(),
    # and since the latter returns assigned centers too, every assigned center
    # was counted twice. Alanine has one chiral center, so the message read 2.
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    taut = _taut("CCC(=O)O", "alanine")

    assert parallel_check_chiral_centers(taut, chirality_facts(contnr)) is None

    message = _logged_message(capsys.readouterr().out)
    assert "total number of chiral centers from 1 to 0" in message
    assert "specified" not in message


def test_chiral_filter_rejects_tautomer_that_drops_a_stereo_assignment() -> None:
    # The artifact tauts_no_elim_chiral was written for: MolVS reports a form
    # that differs from the input only in having dropped a chiral
    # specification. The center survives, so the total count is unchanged (1
    # either way) and a filter comparing only totals would keep it. Keeping it
    # matters because EnumerateChiralMols runs afterwards and expands every
    # unassigned center into both R and S, so the output would carry the
    # enantiomer the input's "@" ruled out.
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    taut = _taut("CC(N)C(=O)O", "alanine")

    assert len(taut.chiral_cntrs_w_unasignd()) == contnr.num_unspecif_chiral_cntrs
    assert parallel_check_chiral_centers(taut, chirality_facts(contnr)) is None


def test_chiral_rejection_names_the_assignment_change(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # A rejection for a changed assignment count must not report it as a change
    # in the total, which is unchanged and would read "from 1 to 1".
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    taut = _taut("CC(N)C(=O)O", "alanine")

    assert parallel_check_chiral_centers(taut, chirality_facts(contnr)) is None

    message = _logged_message(capsys.readouterr().out)
    assert "assigned stereochemistry, from 1 to 0" in message
    assert "total number of chiral centers" not in message


def test_tauts_no_elim_chiral_drops_taut_with_dropped_stereo_assignment() -> None:
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    keep = _taut("C[C@H](N)C(=O)O", "alanine")
    drop = _taut("CC(N)C(=O)O", "alanine")
    result = tauts_no_elim_chiral([contnr], [keep, drop], 1, "serial", None)
    assert result == [keep]


def test_filters_do_not_ship_the_container_to_each_job(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: both filters run once per candidate tautomer and used to pass
    # the whole MolContainer to each job, so every variant of a compound was
    # pickled once per tautomer of it. The comparison needs only a few counts
    # and the SMILES the rejection message quotes.
    payloads: list[tuple[object, ...]] = []

    def capture(params: tuple, num_procs: int, func) -> list:
        """Record the job payloads, then dispatch them as the real path does.

        Args:
            params: The per-job argument tuples the step built.
            num_procs: Processor count, ignored here.
            func: The worker function the step wants applied.

        Returns:
            One result per job, in job order.
        """
        payloads.extend(params)
        return [func(*job) for job in params]

    monkeypatch.setattr(MakeTautomers.Parallelizer, "MultiThreading", capture)

    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {})
    taut = _taut("C[C@H](N)C(=O)O", "alanine")

    assert tauts_no_break_arom_rngs([contnr], [taut], 1, "serial", None) == [taut]
    assert tauts_no_elim_chiral([contnr], [taut], 1, "serial", None) == [taut]

    assert payloads
    for job in payloads:
        assert not any(isinstance(arg, MolContainer) for arg in job)


def test_make_tauts_does_not_ship_the_container_to_each_job(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Same regression in the generation step, which sent the container plus an
    # index into it rather than the molecule, as every other enumeration step
    # does.
    payloads: list[tuple[object, ...]] = []

    def capture(params: tuple, num_procs: int, func) -> list:
        """Record the job payloads, then dispatch them as the real path does.

        Args:
            params: The per-job argument tuples the step built.
            num_procs: Processor count, ignored here.
            func: The worker function the step wants applied.

        Returns:
            One result per job, in job order.
        """
        payloads.extend(params)
        return [func(*job) for job in params]

    monkeypatch.setattr(MakeTautomers.Parallelizer, "MultiThreading", capture)

    contnr = MolContainer("CC(=O)CC", "butanone", 0, {})
    contnr.add_smiles("CC(=O)CC")

    MakeTautomers.make_tauts([contnr], 5, 1, 1, "serial", False, None)

    assert payloads
    for job in payloads:
        assert not any(isinstance(arg, MolContainer) for arg in job)
