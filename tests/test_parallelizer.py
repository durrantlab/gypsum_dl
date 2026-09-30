"""Unit tests for the non-MPI parallelization paths.

MPI mode cannot be exercised without an mpirun launcher, so these tests pin
down the serial and multiprocessing behavior plus the guard rails that fire
when MPI is requested but unavailable.
"""

import multiprocessing
import os
import queue
import random
import sys
import time
import types

import numpy
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


def draw(_ignored: int) -> tuple[float, float]:
    """Draw one number from each generator the worker reseeds.

    Both generators matter: numpy is the one a forked child inherits unchanged,
    and the random module is seeded alongside it so a worker's two streams are
    derived the same way the parent's are. Defined at module scope so worker
    processes can unpickle it.

    Args:
        _ignored: Unused; present so the job takes one argument.

    Returns:
        The number drawn from the random module and the one drawn from numpy.
    """
    return random.random(), float(numpy.random.random())


def numpy_draw_after_delay(_ignored: int) -> float:
    """Draw from the process's numpy legacy global stream after a short pause.

    numpy is the generator that a forked child actually inherits unchanged
    (CPython reseeds the random module through os.register_at_fork, numpy
    installs no such hook), so numpy is what an independence check has to
    sample. The pause keeps each worker busy long enough that the remaining
    tasks go to idle workers rather than piling onto whichever one drained the
    queue first, which is what makes the per-worker streams observable.
    Defined at module scope so worker processes can unpickle it.

    Args:
        _ignored: Unused; present so the job takes one argument.

    Returns:
        The drawn number.
    """
    time.sleep(0.2)
    return float(numpy.random.random())


class _RootComm:
    """Single-rank stand-in for an mpi4py communicator, as the root sees it.

    ParallelMPI.run cannot be exercised under a real launcher from a plain test
    run, and the behavior under test (what happens to a job that raises) is
    entirely on the root's own chunk.
    """

    def Get_rank(self) -> int:
        """Report this rank's index.

        Returns:
            Zero; the stand-in is always the root.
        """
        return 0

    def Get_size(self) -> int:
        """Report the communicator size.

        Returns:
            One, so the root keeps every job.
        """
        return 1

    def bcast(self, obj: object, root: int = 0) -> object:
        """Hand the broadcast object straight back.

        Args:
            obj: The object being broadcast.
            root: Ignored; present to match the mpi4py signature.

        Returns:
            The object unchanged.
        """
        return obj

    def scatter(self, obj: list[object], root: int = 0) -> object:
        """Return the chunk destined for rank zero.

        Args:
            obj: The per-rank chunks.
            root: Ignored; present to match the mpi4py signature.

        Returns:
            The first chunk.
        """
        return obj[0]

    def gather(self, obj: object, root: int = 0) -> list[object]:
        """Collect this rank's results.

        Args:
            obj: This rank's result chunk.
            root: Ignored; present to match the mpi4py signature.

        Returns:
            A one-element list of rank results.
        """
        return [obj]


class _WorkerComm:
    """Stand-in communicator that feeds one chunk to ParallelMPI._worker.

    The worker loop is an infinite receive loop, so the stand-in also has to be
    able to end it: the second broadcast is the kill signal.
    """

    def __init__(
        self,
        func: object,
        chunk: list[list[object]],
        seed_chunk: list[list[int]],
    ) -> None:
        """Record what the worker should receive.

        Args:
            func: The job function delivered by the first broadcast.
            chunk: The argument chunk delivered by the first scatter.
            seed_chunk: The per-job seeds delivered by the second scatter.
        """
        self._bcasts: list[object] = [func, None]
        self._scatters: list[object] = [chunk, seed_chunk]
        self.gathered: list[object] = []

    def bcast(self, obj: object, root: int = 0) -> object:
        """Deliver the next scripted broadcast.

        Args:
            obj: Ignored; the worker always passes None.
            root: Ignored; present to match the mpi4py signature.

        Returns:
            The job function, then None to stop the loop.
        """
        return self._bcasts.pop(0)

    def scatter(self, obj: object, root: int = 0) -> object:
        """Deliver this worker's next scattered chunk.

        Args:
            obj: Ignored; the worker always passes an empty list.
            root: Ignored; present to match the mpi4py signature.

        Returns:
            The argument chunk, then the matching seed chunk.
        """
        return self._scatters.pop(0)

    def gather(self, obj: object, root: int = 0) -> None:
        """Record what the worker tried to send back.

        Args:
            obj: This worker's result chunk.
            root: Ignored; present to match the mpi4py signature.

        Returns:
            None, as the non-root return value of an mpi4py gather.
        """
        self.gathered.append(obj)
        return None


def _make_parallel_mpi(
    monkeypatch: pytest.MonkeyPatch, comm: object
) -> parallelizer.ParallelMPI:
    """Build a ParallelMPI whose communicator is a stand-in.

    The constructor reads mpi4py.MPI.COMM_WORLD out of the module namespace,
    and mpi4py is not importable in a plain test run, so the module reference
    is what has to be replaced.

    Args:
        monkeypatch: Fixture used to install the stand-in module.
        comm: The stand-in communicator.

    Returns:
        A ParallelMPI bound to `comm`.
    """
    stub = types.ModuleType("mpi4py")
    stub_mpi = types.ModuleType("mpi4py.MPI")
    stub_mpi.COMM_WORLD = comm
    stub.MPI = stub_mpi
    monkeypatch.setattr(parallelizer, "mpi4py", stub, raising=False)
    return parallelizer.ParallelMPI()


def _drain_worker(seeds: list[int]) -> list[tuple[float, float]]:
    """Run parallelizer.worker in this process against a prefilled queue.

    Running the worker body in-process is the only way to observe its seeding
    deterministically: what a forked child's generator does is otherwise
    visible only through the results of a race between workers. Plain
    queue.Queue is enough, since the worker only ever calls get and put.

    Args:
        seeds: One seed per job, carried by the job the way MultiThreading
            attaches them.

    Returns:
        The (random, numpy) draw pairs, in job order.
    """
    task_queue: queue.Queue = queue.Queue()
    done_queue: queue.Queue = queue.Queue()
    for index, seed in enumerate(seeds):
        task_queue.put((index, (draw, (index,), seed)))
    task_queue.put("STOP")

    parallelizer.worker(task_queue, done_queue)

    results = [done_queue.get() for _ in range(len(seeds))]
    results.sort(key=lambda pair: pair[0])
    return [pair[1] for pair in results]


def test_run_one_passes_arguments_through() -> None:
    assert parallelizer.run_one(add_one, (1,)) == 2


def test_run_one_returns_none_on_exception() -> None:
    assert parallelizer.run_one(add_one_unless_two, (2,)) is None


def test_run_one_names_the_failed_function(capsys: pytest.CaptureFixture[str]) -> None:
    parallelizer.run_one(add_one_unless_two, (2,))
    captured = capsys.readouterr().out
    assert "add_one_unless_two" in captured
    assert "simulated job failure" in captured


def test_run_one_names_the_molecule_that_failed(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Regression: the report named the function and nothing else, so a
    # production log recorded that a step had failed without recording which
    # compound it failed on.
    class _Subject:
        orig_smi = "CC(=O)C"
        name = "acetone"

    def raiser(subject: object) -> None:
        """Fail the way a real job does, with the subject as its argument.

        Args:
            subject: The molecule the job was working on.

        Raises:
            RuntimeError: Always.
        """
        raise RuntimeError("simulated job failure")

    assert parallelizer.run_one(raiser, (_Subject(),)) is None

    captured = capsys.readouterr().out
    assert "CC(=O)C (acetone)" in captured
    assert "raiser" in captured


def test_describe_job_subject_reads_a_mapping() -> None:
    # The tautomer step passes its container-level fields as a mapping rather
    # than shipping the container to every job.
    props = {"orig_smi": "CCO", "name": "ethanol"}
    assert parallelizer.describe_job_subject((props,)) == "CCO (ethanol)"


def test_describe_job_subject_tolerates_unidentifiable_arguments() -> None:
    # This runs inside an exception handler, and a half-built object is a
    # common reason for a job to have raised, so it must not raise itself and
    # replace the traceback the user needs.
    assert parallelizer.describe_job_subject(()) == ""
    assert parallelizer.describe_job_subject((7,)) == ""
    assert parallelizer.describe_job_subject(({"name": "ethanol"},)) == "ethanol"


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


@pytest.mark.parametrize("version", ["4.0.1", "4.0", "4.1.0rc1", "2", "10.0"])
def test_mpi4py_version_supported_accepts_usable_releases(version: str) -> None:
    # Regression: the version was parsed by running int() over every
    # dot-separated component and indexing the minor unconditionally, so a
    # suffixed release ("4.1.0rc1") raised ValueError and a single-component
    # version ("2") raised IndexError. An unusual version string is not
    # evidence of an old mpi4py.
    assert parallelizer.mpi4py_version_supported(version) is True


@pytest.mark.parametrize("version", ["3.1.4", "2.1.0", "2.1", "2.0.1", "1.3.1"])
def test_mpi4py_version_supported_rejects_old_releases(version: str) -> None:
    # Regression: the gate was set at (2, 1), the mpi4py release that added
    # the "-m mpi4py" launch flag, while pyproject.toml and pixi.toml both
    # require mpi4py>=4.0.1. It could therefore never fire, and its message
    # advertised support for releases nothing tests against.
    assert parallelizer.mpi4py_version_supported(version) is False


def test_mpi4py_version_gate_matches_the_declared_dependency() -> None:
    # The gate and the packaging pins have to move together, and the message
    # has to quote the gate rather than a version of its own.
    assert parallelizer.MIN_MPI4PY_VERSION == (4, 0)
    assert "4.0 or higher" in parallelizer.MPI4PY_VERSION_MSG


@pytest.mark.parametrize(
    ("version", "expected"),
    [("4.1.0rc1", True), ("2", True), ("4.0.1", True), ("3.1.4", False)],
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


def test_parallel_mpi_run_drops_a_failed_job(monkeypatch: pytest.MonkeyPatch) -> None:
    # Regression: the root's own chunk called func(*arg) bare, so one raising
    # molecule propagated out of run() before COMM.gather, leaving every other
    # rank blocked in gather until the scheduler killed the job. This is the mpi
    # twin of test_multithreading_serial_path_drops_failed_job.
    par = _make_parallel_mpi(monkeypatch, _RootComm())

    assert par.run(add_one_unless_two, [[1], [2], [3]]) == [2, None, 4]


def test_parallel_mpi_worker_drops_a_failed_job(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: the worker's chunk had the same bare call, so a molecule that
    # raised on a non-root rank skipped that rank's gather and hung the run.
    comm = _WorkerComm(add_one_unless_two, [[1], [2], [3]], [[11], [12], [13]])
    par = _make_parallel_mpi(monkeypatch, comm)

    # The worker loop only ends on the kill signal, which exits the process.
    with pytest.raises(SystemExit):
        par._worker()

    assert comm.gathered == [[2, None, 4]]


def test_parallel_mpi_run_seeds_jobs_like_the_other_modes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The root draws one seed per job and scatters the seeds alongside the
    # arguments, so a seeded mpi run has to land on the same numbers a seeded
    # serial run does.
    inputs = [(i,) for i in range(4)]
    parallelizer.seed_generators(99)
    serial = parallelizer.MultiThreading(inputs, 1, draw)

    par = _make_parallel_mpi(monkeypatch, _RootComm())
    parallelizer.seed_generators(99)

    assert par.run(draw, [list(item) for item in inputs]) == serial


def test_parallel_mpi_worker_applies_the_scattered_seeds(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # A non-root rank has no way to know a job's position in the input, so its
    # seeds have to arrive with the work. Without them the rank would draw from
    # whatever state it was left in when it parked in the worker loop.
    parallelizer.seed_generators(99)
    serial = parallelizer.MultiThreading([(i,) for i in range(3)], 1, draw)

    parallelizer.seed_generators(99)
    seed_chunk = [[seed] for seed in parallelizer.draw_job_seeds(3)]
    comm = _WorkerComm(draw, [[0], [1], [2]], seed_chunk)
    par = _make_parallel_mpi(monkeypatch, comm)

    with pytest.raises(SystemExit):
        par._worker()

    assert comm.gathered == [serial]


def test_parallel_mpi_worker_skips_the_padding_split_hands_it(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The worker's dedicated filler branch compared an argument list to
    # Empty_obj and could never fire, so padding was skipped only by an
    # incidental check further down. Taking the chunk from _split, rather than
    # spelling out its shape here, keeps this test tied to the padding the
    # root actually sends.
    par = _make_parallel_mpi(monkeypatch, _RootComm())
    arg_chunks = par._split([[1]], 3)
    seed_chunks = par._split([[11]], 3)

    comm = _WorkerComm(add_one_unless_two, arg_chunks[-1], seed_chunks[-1])
    par.COMM = comm

    with pytest.raises(SystemExit):
        par._worker()

    assert comm.gathered == [[]]


def test_parallel_mpi_failure_handling_matches_the_other_modes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    inputs = [(1,), (2,), (3,)]
    serial_par = parallelizer.Parallelizer("serial", 1)
    serial = serial_par.run(inputs, add_one_unless_two)
    serial_par.end()

    mpi_par = _make_parallel_mpi(monkeypatch, _RootComm())

    assert mpi_par.run(add_one_unless_two, [list(i) for i in inputs]) == serial


def test_worker_applies_the_seed_carried_by_each_job() -> None:
    # A job's stream has to follow from its own seed and nothing else: not from
    # the state the worker started with, and not from how many jobs it ran
    # before this one. Two jobs given the same seed therefore draw the same
    # numbers, whatever state the process was in beforehand.
    random.seed(1)
    numpy.random.seed(1)
    first = _drain_worker([7, 8, 7])

    random.seed(12345)
    numpy.random.seed(12345)
    again = _drain_worker([7, 8, 7])

    assert first == again
    assert first[0] == first[2]
    assert first[0] != first[1]


def test_run_one_without_a_seed_draws_from_the_state_in_place() -> None:
    # Each pipeline step calls run_one directly when it was handed no
    # parallelizer object, and those calls have no seed of their own to
    # install, so the unseeded path has to keep drawing from the state the
    # caller set up.
    random.seed(3)
    numpy.random.seed(3)
    expected = draw(0)

    random.seed(3)
    numpy.random.seed(3)

    assert parallelizer.run_one(draw, (0,)) == expected


def test_run_one_restores_the_callers_generator_state() -> None:
    # A seeded job that runs in the dispatching process (serial mode, the mpi
    # root's own chunk) has to leave the caller's streams where a job run in a
    # child process would leave them. Otherwise the seeds drawn for the next
    # step would depend on the job manager.
    random.seed(3)
    numpy.random.seed(3)
    expected = draw(0)

    random.seed(3)
    numpy.random.seed(3)
    parallelizer.run_one(draw, (0,), 424242)

    assert draw(0) == expected


def test_draw_job_seeds_is_reproducible_and_distinct() -> None:
    parallelizer.seed_generators(5)
    first = parallelizer.draw_job_seeds(6)

    parallelizer.seed_generators(5)

    assert parallelizer.draw_job_seeds(6) == first
    assert len(set(first)) == 6


def test_seeded_results_do_not_depend_on_the_number_of_workers() -> None:
    # Regression: the seeds were drawn one per worker rather than one per job,
    # which cannot make a run reproducible. The task queue is drained by
    # whichever worker is free, so a molecule's stream depended on who picked
    # it up, and the run varied between invocations while looking as though a
    # seed had pinned it down.
    inputs = [(i,) for i in range(8)]

    runs = []
    for num_procs in (1, 2, 3, 8):
        parallelizer.seed_generators(1234)
        runs.append(parallelizer.MultiThreading(inputs, num_procs, draw))

    assert runs[0] == runs[1] == runs[2] == runs[3]

    # Reproducible, but still a separate stream per job.
    assert len(set(runs[0])) == len(inputs)


def test_seeded_results_change_with_the_seed() -> None:
    inputs = [(i,) for i in range(4)]

    parallelizer.seed_generators(1234)
    first = parallelizer.MultiThreading(inputs, 2, draw)

    parallelizer.seed_generators(4321)
    second = parallelizer.MultiThreading(inputs, 2, draw)

    assert first != second


def test_worker_processes_draw_independent_numpy_streams() -> None:
    # Regression: start_processes forks the workers after the parent has seeded
    # itself, and numpy installs no os.register_at_fork hook, so every child
    # held a byte-identical copy of the numpy legacy global state and nothing
    # reseeded it. The ring-conformation clustering picks its initial centroids
    # from that generator (scipy's kmeans2 with minit="points"), so every
    # worker clustered against the same draw. Nothing about the output looks
    # wrong; the conformer selection is just correlated across workers.
    #
    # Note that the random module needs no such fix: CPython reseeds the global
    # instance in the forked child, so a test written against random.random()
    # passes either way and proves nothing.
    #
    # The seed now travels with the job rather than with the worker, so the
    # streams are distinct per job; the delay in the job is left in place so
    # that the tasks really do spread across the workers.
    draws = parallelizer.MultiThreading(
        [(i,) for i in range(4)], 4, numpy_draw_after_delay
    )

    assert len(set(draws)) == 4


@pytest.mark.parametrize(
    "orig_argv",
    [
        ["python", "-m", "mpi4py", "-m", "gypsum_dl"],
        ["python", "-m", "mpi4py", "/env/bin/gypsum-dl", "-j", "p.json"],
        ["python", "-u", "-X", "dev", "-W", "ignore", "-m", "mpi4py", "-m", "x"],
        ["python", "-um", "mpi4py", "-m", "gypsum_dl"],
        ["python", "-mmpi4py", "-m", "gypsum_dl"],
        ["python", "-Wignore", "-m", "mpi4py", "-m", "gypsum_dl"],
    ],
)
def test_mpi4py_launch_flag_present_accepts_mpi4py_launches(
    monkeypatch: pytest.MonkeyPatch, orig_argv: list[str]
) -> None:
    monkeypatch.setattr(sys, "orig_argv", orig_argv)
    assert parallelizer.mpi4py_launch_flag_present() is True


@pytest.mark.parametrize(
    "orig_argv",
    [
        ["python", "-m", "gypsum_dl", "--job_manager", "mpi"],
        ["/env/bin/python", "/env/bin/gypsum-dl", "--job_manager", "mpi"],
        ["python", "script.py", "-m", "mpi4py"],
        ["python", "-c", "import x", "-m", "mpi4py"],
        ["python", "-W", "-m", "mpi4py"],
        ["python", "-m"],
        [],
    ],
)
def test_mpi4py_launch_flag_present_rejects_other_launches(
    monkeypatch: pytest.MonkeyPatch, orig_argv: list[str]
) -> None:
    # Regression: the flag was inferred from runpy being in sys.modules, which
    # any other route that loads runpy also satisfies. "python -m gypsum_dl"
    # with --job_manager mpi passed the check without mpi4py's exception
    # handling in place, so a failing rank hung the job instead of aborting it.
    monkeypatch.setattr(sys, "orig_argv", orig_argv)
    monkeypatch.setitem(sys.modules, "runpy", types.ModuleType("runpy"))
    assert parallelizer.mpi4py_launch_flag_present() is False
