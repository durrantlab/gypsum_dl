"""Unit tests for the non-MPI parallelization paths.

MPI mode cannot be exercised without an mpirun launcher, so these tests pin
down the serial and multiprocessing behavior plus the guard rails that fire
when MPI is requested but unavailable.
"""

import multiprocessing

import pytest

from gypsum_dl import parallelizer


def add_one(value: int) -> int:
    """Increment a value.

    Defined at module scope so worker processes can unpickle it.

    Args:
        value: Number to increment.

    Returns:
        The incremented value.
    """
    return value + 1


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