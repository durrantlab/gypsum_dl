"""Unit tests for the non-MPI parallelization paths.

MPI mode cannot be exercised without an mpirun launcher, so these tests pin
down the serial and multiprocessing behavior plus the guard rails that fire
when MPI is requested but unavailable.
"""

import multiprocessing
import os
import sys
import time
import types

import pytest

from gypsum_dl import parallelizer


def _stub_mpi4py(monkeypatch: pytest.MonkeyPatch, version: str) -> None:
    """Install a stand-in mpi4py reporting a chosen version.

    check_mpi_available imports both mpi4py and its MPI sublibrary, neither of
    which is available (or launchable) in a plain test run, so an arbitrary
    version string can only be exercised through stand-ins.

    Args:
        monkeypatch: Fixture used to install the stand-ins.
        version: The version string the stand-in should report.
    """
    stub = types.ModuleType("mpi4py")
    stub.__version__ = version
    stub_mpi = types.ModuleType("mpi4py.MPI")
    stub.MPI = stub_mpi
    monkeypatch.setitem(sys.modules, "mpi4py", stub)
    monkeypatch.setitem(sys.modules, "mpi4py.MPI", stub_mpi)


def add_one(value: int) -> int:
    """Increment a value.

    Defined at module scope so worker processes can unpickle it.

    Args:
        value: Number to increment.

    Returns:
        The incremented value.
    """
    return value + 1


def add_one_unless_two(value: int) -> int:
    """Increment a value, but blow up on one specific input.

    Stands in for a job that hits an edge case in only one of its inputs.
    Defined at module scope so worker processes can unpickle it.

    Args:
        value: Number to increment.

    Returns:
        The incremented value.

    Raises:
        ValueError: If `value` is 2.
    """
    if value == 2:
        raise ValueError("simulated job failure")
    return value + 1


def kill_worker_on_two(value: int) -> int:
    """Increment a value, but kill the calling process on one specific input.

    Stands in for a job whose worker dies without reporting a result, the way
    an OOM kill or a native fault in RDKit does. os._exit is used instead of
    raising so that no exception handler can turn it back into a None result.
    Defined at module scope so worker processes can unpickle it.

    Args:
        value: Number to increment.

    Returns:
        The incremented value, for every input but 2.
    """
    if value == 2:
        os._exit(1)
    return value + 1


def test_run_one_passes_arguments_through() -> None:
    assert parallelizer.run_one(add_one, (1,)) == 2


def test_run_one_returns_none_on_exception() -> None:
    assert parallelizer.run_one(add_one_unless_two, (2,)) is None


def test_run_one_names_the_failed_function(capsys: pytest.CaptureFixture[str]) -> None:
    parallelizer.run_one(add_one_unless_two, (2,))
    captured = capsys.readouterr().out
    assert "add_one_unless_two" in captured
    assert "simulated job failure" in captured


def test_flatten_list_handles_none() -> None:
    assert parallelizer.flatten_list(None) == []


def test_flatten_list_nested() -> None:
    assert parallelizer.flatten_list([[1, 2], [3]]) == [1, 2, 3]


def test_flatten_list_already_flat() -> None:
    assert parallelizer.flatten_list([1, 2]) == [1, 2]


def test_flatten_list_drops_none_among_sublists() -> None:
    # A None worker result mixed with lists must be dropped, not returned
    # unflattened. Regression for flatten_list becoming a no-op.
    assert parallelizer.flatten_list([[1, 2], None, [3]]) == [1, 2, 3]


def test_flatten_list_drops_none_within_sublist() -> None:
    assert parallelizer.flatten_list([[1, None, 2], [3]]) == [1, 2, 3]


def test_strip_none_removes_nones() -> None:
    assert parallelizer.strip_none([1, None, 2]) == [1, 2]


def test_strip_none_handles_none_input() -> None:
    assert parallelizer.strip_none(None) == []


def test_count_processors_caps_at_number_of_inputs() -> None:
    assert parallelizer.count_processors(2, 8) == 2


def test_count_processors_uses_all_cpus_when_non_positive() -> None:
    assert parallelizer.count_processors(10**6, 0) == multiprocessing.cpu_count()


def test_check_and_format_inputs_converts_lists_to_tuples() -> None:
    assert parallelizer.check_and_format_inputs_to_list_of_tuples([[1], [2]]) == [
        (1,),
        (2,),
    ]


def test_check_and_format_inputs_passes_tuples_through() -> None:
    args = [(1,), (2,)]
    assert parallelizer.check_and_format_inputs_to_list_of_tuples(args) is args


def test_check_and_format_inputs_rejects_non_sequence() -> None:
    with pytest.raises(Exception, match="list of tuples"):
        parallelizer.check_and_format_inputs_to_list_of_tuples("nope")


def test_check_and_format_inputs_rejects_mixed_types() -> None:
    with pytest.raises(Exception, match="same type"):
        parallelizer.check_and_format_inputs_to_list_of_tuples([(1,), [2]])


def test_multithreading_empty_input() -> None:
    assert parallelizer.MultiThreading([], 4, add_one) == []


def test_multithreading_serial_path() -> None:
    assert parallelizer.MultiThreading([(1,), (2,)], 1, add_one) == [2, 3]


def test_multithreading_preserves_input_order_across_processes() -> None:
    inputs = [(i,) for i in range(4)]
    assert parallelizer.MultiThreading(inputs, 2, add_one) == [1, 2, 3, 4]


def test_multithreading_serial_path_drops_failed_job() -> None:
    # Regression: the single-processor branch used to call the job bare, so a
    # raising input aborted the whole run instead of yielding None.
    inputs = [(1,), (2,), (3,)]
    assert parallelizer.MultiThreading(inputs, 1, add_one_unless_two) == [2, None, 4]


def test_multithreading_raises_when_a_worker_dies(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: the collection loop blocked on exactly len(inputs)
    # done_queue.get() calls with no timeout and no liveness check, so a worker
    # that died without reporting its result left the parent waiting forever.
    # Shorten the poll interval so the test does not sit through the
    # production-sized wait.
    monkeypatch.setattr(parallelizer, "WORKER_POLL_TIMEOUT", 1.0)

    inputs = [(i,) for i in range(4)]
    started = time.monotonic()
    with pytest.raises(Exception, match="died"):
        parallelizer.MultiThreading(inputs, 3, kill_worker_on_two)

    # The point of the fix is that the failure is bounded rather than a hang.
    assert time.monotonic() - started < 60.0


def test_multithreading_leaves_no_orphan_children() -> None:
    # Regression: worker processes were started but never joined, so a
    # completed run left non-daemon children behind.
    parallelizer.MultiThreading([(i,) for i in range(4)], 2, add_one)
    assert multiprocessing.active_children() == []


def test_multithreading_failure_handling_matches_across_procs() -> None:
    inputs = [(1,), (2,), (3,)]
    serial = parallelizer.MultiThreading(inputs, 1, add_one_unless_two)
    parallel = parallelizer.MultiThreading(inputs, 2, add_one_unless_two)
    assert serial == parallel


def test_parallelizer_failure_handling_matches_across_modes() -> None:
    inputs = [(1,), (2,), (3,)]

    serial_par = parallelizer.Parallelizer("serial", 1)
    serial = serial_par.run(inputs, add_one_unless_two)
    serial_par.end()

    mp_par = parallelizer.Parallelizer("multiprocessing", 2, True)
    parallel = mp_par.run(inputs, add_one_unless_two)
    mp_par.end()

    assert serial == parallel == [2, None, 4]


def test_parallelizer_serial_mode_runs_jobs() -> None:
    par = parallelizer.Parallelizer("serial", 4)
    assert par.return_mode() == "serial"
    assert par.return_node() == 1
    assert par.run([(1,), (2,)], add_one) == [2, 3]
    par.end()


def test_parallelizer_multiprocessing_mode_runs_jobs() -> None:
    par = parallelizer.Parallelizer("multiprocessing", 1, True)
    assert par.return_mode() == "multiprocessing"
    assert par.run([(1,), (2,)], add_one) == [2, 3]
    par.end()


def test_parallelizer_defaults_num_procs_to_cpu_count() -> None:
    par = parallelizer.Parallelizer("multiprocessing", None, True)
    assert par.return_node() == multiprocessing.cpu_count()


def test_parallelizer_negative_num_procs_falls_back_to_cpu_count() -> None:
    par = parallelizer.Parallelizer("multiprocessing", -1, True)
    assert par.return_node() == multiprocessing.cpu_count()


def test_parallelizer_rejects_mpi_when_unavailable() -> None:
    with pytest.raises(Exception, match="mpi4py"):
        parallelizer.Parallelizer("mpi", 1, True)


def test_parallelizer_run_cannot_switch_to_mpi() -> None:
    par = parallelizer.Parallelizer("multiprocessing", 1, True)
    with pytest.raises(Exception, match="non-mpi to mpi"):
        par.run([(1,)], add_one, mode="mpi")


def test_parallelizer_run_rejects_unknown_mode() -> None:
    par = parallelizer.Parallelizer("multiprocessing", 1, True)
    with pytest.raises(Exception, match="doesn't match"):
        par.run([(1,)], add_one, mode="bogus")


def test_parallelizer_run_logs_the_mode_change_in_the_right_direction(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Regression: the message was formatted as "from {mode} to {self.mode}",
    # naming the destination as the source. run() dispatches on `mode`, so the
    # override direction is self.mode -> mode.
    par = parallelizer.Parallelizer("multiprocessing", 1, True)

    par.run([(1,)], add_one, mode="serial")
    par.end()

    assert "changing mode from multiprocessing to serial" in capsys.readouterr().out


def test_parallelizer_run_rejects_num_procs_override_in_serial_mode() -> None:
    par = parallelizer.Parallelizer("serial", 1)
    with pytest.raises(Exception, match="Can't override num_procs"):
        par.run([(1,)], add_one, num_procs=4, mode="serial")


def test_parallelizer_compute_nodes_per_mode() -> None:
    par = parallelizer.Parallelizer("serial", 1)
    assert par.compute_nodes("serial") == 1
    assert par.compute_nodes("multiprocessing") == multiprocessing.cpu_count()
    with pytest.raises(Exception, match="mpi4py"):
        par.compute_nodes("mpi")


def test_parallelizer_start_returns_none_for_non_mpi() -> None:
    par = parallelizer.Parallelizer("serial", 1)
    assert par.start("multiprocessing") is None
    with pytest.raises(Exception, match="mpi4py"):
        par.start("mpi")


def test_parallelizer_end_rejects_mpi_when_unavailable() -> None:
    par = parallelizer.Parallelizer("serial", 1)
    with pytest.raises(Exception, match="mpi4py"):
        par.end("mpi")


def test_parallelizer_pick_mode_method_survives_init() -> None:
    # Regression: __init__ stored the result of pick_mode() back onto the
    # attribute `self.pick_mode`, clobbering the bound method with a string.
    # Any later self.pick_mode() call would then raise TypeError. The picked
    # mode now lives on `picked_mode`, leaving the method callable.
    par = parallelizer.Parallelizer(None, 1, True)
    assert callable(par.pick_mode)
    assert par.picked_mode == "multiprocessing"
    assert par.pick_mode() == par.picked_mode
    par.end()


@pytest.mark.parametrize("version", ["2.1.0", "2.1", "3.1.4", "4.1.0rc1", "2", "10.0"])
def test_mpi4py_version_supported_accepts_usable_releases(version: str) -> None:
    # Regression: the version was parsed by running int() over every
    # dot-separated component and indexing the minor unconditionally, so a
    # suffixed release ("4.1.0rc1") raised ValueError and a single-component
    # version ("2") raised IndexError. An unusual version string is not
    # evidence of an old mpi4py.
    assert parallelizer.mpi4py_version_supported(version) is True


@pytest.mark.parametrize("version", ["2.0.1", "1.3.1"])
def test_mpi4py_version_supported_rejects_old_releases(version: str) -> None:
    assert parallelizer.mpi4py_version_supported(version) is False


@pytest.mark.parametrize(
    ("version", "expected"),
    [("4.1.0rc1", True), ("2", True), ("3.1.4", True), ("2.0.1", False)],
)
def test_check_mpi_available_follows_the_shared_version_gate(
    monkeypatch: pytest.MonkeyPatch, version: str, expected: bool
) -> None:
    # Regression: this check had its own copy of the version parse, so the
    # version strings that start handles fine raised in here instead. The
    # caller catches every exception and returns False, which silently demoted
    # a perfectly good mpi4py to multiprocessing.
    _stub_mpi4py(monkeypatch, version)
    par = parallelizer.Parallelizer("serial", 1)

    assert par.check_mpi_available() is expected


def test_mpi4py_launch_flag_present_reads_sys_modules(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setitem(sys.modules, "runpy", types.ModuleType("runpy"))
    assert parallelizer.mpi4py_launch_flag_present() is True

    monkeypatch.delitem(sys.modules, "runpy")
    assert parallelizer.mpi4py_launch_flag_present() is False
