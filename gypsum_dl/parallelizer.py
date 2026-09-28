"""
Parallelizer.py

Abstract parallel computation utility.

The "parallelizer" object exposes a simple map interface that takes a function
and a list of arguments and returns the result of applying the function to each
argument. Internally, the parallelizer class can determine what parallel
capabilities are present on a system and automatically pick between "mpi",
"multiprocessing" or "serial" in order to speed up the map operation. This
approach simplifies development and allows the same program to run on a laptop
or a high-performance computer cluster, utilizing the full resources of each
system. (Description provided by Harrison Green.)
"""

from typing import Any, TypeVar

import multiprocessing
import queue
import random
import sys
import traceback
from collections.abc import Callable, Sequence

import numpy

_JobResult = TypeVar("_JobResult")

WORKER_POLL_TIMEOUT: float = 60.0
"""Seconds to wait for a worker result before checking whether any worker is
still alive. Long jobs are fine: the wait simply repeats while at least one
worker is running."""

WORKER_EXIT_GRACE: float = 1.0
"""Seconds to keep draining results after the last worker exits, so a result
that was in flight when the worker finished is not mistaken for a lost one."""

try:
    import mpi4py

    MPI_installed = True
except Exception:
    MPI_installed = False

MIN_MPI4PY_VERSION: tuple[int, int] = (4, 0)
"""Oldest mpi4py gypsum-dl runs against. The "-m mpi4py" launch flag the mpi
job manager depends on arrived in 2.1.0, but the declared dependency in
pyproject.toml and pixi.toml is mpi4py>=4.0.1, so a gate set at the old
feature floor could never fire and advertised support for releases nobody
tests against. Keep this in step with those pins."""

MPI_LAUNCH_FLAG_MSG: str = (
    "\nTo run in mpi mode you must run with -m flag. ie) mpirun -n $NTASKS python -m mpi4py run_gypsum_dl.py\n"
)

MPI4PY_MISSING_MSG: str = (
    "\nmpi4py not installed but --job_manager is set to mpi. \n Either install mpi4py or switch job_manager to multiprocessing or serial.\n"
)

MPI4PY_VERSION_MSG: str = (
    "\nmpi4py version "
    + ".".join(str(part) for part in MIN_MPI4PY_VERSION)
    + " or higher is required. Use the 'python -m mpi4py' flag to run in mpi mode.\nPlease update mpi4py to a newer version, or switch job_manager to multiprocessing or serial.\n"
)


def mpi4py_launch_flag_present() -> bool:
    """Report whether python was started with the "-m mpi4py" runpy flag.

    mpi4py overrides the way exceptions propagate, so this has to be settled
    before the api is loaded. Both the parameter setup in start and the
    Parallelizer's own probe ask the question, and they have to agree on how it
    is detected.

    Returns:
        True if the flag appears to have been used.
    """

    return "runpy" in sys.modules


def mpi4py_version_supported(version: str) -> bool:
    """Judge an mpi4py version string against the minimum gypsum-dl needs.

    Parses only the two components the comparison needs. Running int() over
    every dot-separated component raised ValueError on a suffixed release
    ("4.1.0rc1"), and indexing the minor unconditionally raised IndexError on a
    single-component version ("2"), so a check meant to produce a friendly
    message instead ended the run with a traceback (or, in the Parallelizer,
    quietly demoted a perfectly good mpi4py to multiprocessing).

    Args:
        version: The version string reported by mpi4py.

    Returns:
        True if the version is recent enough, or if it cannot be parsed: an
            unusual version string is not evidence of an old mpi4py.
    """

    try:
        major, minor = (int(x) for x in version.split(".")[:2])
    except ValueError:
        return True

    return (major, minor) >= MIN_MPI4PY_VERSION


class Parallelizer(object):
    """
    Abstract parallelization class
    """

    def __init__(
        self,
        mode: str | None = None,
        num_procs: int | None = None,
        flag_for_low_level: bool = False,
    ):
        """
        This will initialize the Parallelizer class and kick off the specific
        classes for multiprocessing and MPI.

        Default num_procs is all the processesors possible.

        Args:
            mode: the multiprocess mode to be used, ie) serial,
                multiprocessing, mpi, or None: if None then we will try to
                pick a possible multiprocessing choice. This should only
                be used for top level coding. It is best practice to specify
                which multiprocessing choice to use. if you have smaller programs
                used by a larger program, with both mpi enabled there will be
                problems, so specify multiprocessing is important.
            num_procs: the number of processors or nodes that will be used.
                If None than we will use all available nodes/processors
                This will be overriden and fixed to a single processor if
                mode==serial.
            flag_for_low_level: this will override mode and number of processors
                and set it to a multiprocess as serial. This is useful because
                a low-level program in mpi mode referenced by a top level program
                in mpi mode will have terrible problems. This means you can't
                mpi-multiprocess inside an mpi-multiprocess.

        Attrs:
            mode:  The mode which will be used for this paralellization.
                This determined by mode and the enviorment default is mpi if the
                enviorment allows, if not mpi then multiprocessing on all
                processors unless stated otherwise.
            parallel_obj: This is the obstantiated object of the class of parallizations.
                ie) if self.mode=='mpi' self.parallel_obj will be an instance of the mpi class
                This style of retained parallel_obj to be used later is important because
                this is the object which controls the work nodes and maintains the mpi universe
                self.parallel_obj will be set to None for simpler parallization methods like serial
            num_processor: the number of processors or nodes that will be used. If None than we
                will use all available nodes/processors This will be overriden and fixed to a
                single processor if mode==serial
        """

        if mode in ["none", "None"]:
            mode = None

        self.HAS_MPI: bool = self.test_import_MPI(mode, flag_for_low_level)
        """
        If true it can import mpi4py; if False it either cant import mpi4py or the
        mode/flag_for_low_level dictates not to check mpi4py import due to issues
        which arrise when mpi is started within a program which is already mpi enabled.
        """

        # Pick the mode
        self.picked_mode = self.pick_mode()

        if mode is None:
            if self.picked_mode == "mpi" and self.HAS_MPI == True:
                # THIS IS TO BE RUN IN MPI
                self.mode = "mpi"
            else:
                self.mode = "multiprocessing"

        elif mode == "mpi":
            if self.HAS_MPI == True:
                # THIS IS EXPLICITILY CHOSEN TO BE RUN IN MPI AND CAN WORK WITH MPI
                self.mode = "mpi"
            else:
                raise Exception("mpi4py package must be available to use mpi mode")

        elif mode == "multiprocessing":
            self.mode = "multiprocessing"

        elif mode in ["Serial", "serial"]:
            self.mode = "serial"

        else:
            # Default setting will be multiprocessing
            self.mode = "multiprocessing"

        # Start MPI MODE if applicable
        self.parallel_obj = self.start(self.mode) if self.mode == "mpi" else None

        if self.mode == "serial":
            self.num_procs = 1

        elif num_procs == None:
            self.num_procs = self.compute_nodes()

        elif num_procs >= 1:
            self.num_procs = num_procs

        else:
            self.num_procs = self.compute_nodes()

    def test_import_MPI(self, mode, flag_for_low_level=False):
        """
        This tests for the ability of importing the MPI sublibrary from mpi4py.

        This import is problematic when run inside a program which was already mpi parallelized (ie a program run inside an mpi program)
            - for some reason from mpi4py import MPI is problematic in this sub-program structuring.
        To prevent these errors we do a quick check outside the class with a Try statement to import mpi4py
            - if it can't do the import mpi4py than the API isn't installed and we can't run MPI, so we won't even attempt from mpi4py import MPI

        it then checks if the mode has been already establish or if there is a low level flag.

        If the user explicitly or implicitly asks for mpi  (ie mode=None or mode="mpi") without flags and mpi4py is installed, then we will
            run the from mpi4py import MPI check. if it passes then we will return a True and run mpi mode; if not we return False and run multiprocess

        Inputs:
        :param str mode: the multiprocess mode to be used, ie) serial, multiprocessing, mpi, or None:
                            if None then we will try to pick a possible multiprocessing choice. This should only be used for
                            top level coding. It is best practice to specify which multiprocessing choice to use.
                            if you have smaller programs used by a larger program, with both mpi enabled there will be problems, so specify multiprocessing is important.
        :param int num_procs:   the number of processors or nodes that will be used. If None than we will use all available nodes/processors
                                        This will be overriden and fixed to a single processor if mode==serial
        :param bol flag_for_low_level: this will override mode and number of processors and set it to a multiprocess as serial. This is useful because
                                a low-level program in mpi mode referenced by a top level program in mpi mode will have terrible problems. This means you can't mpi-multiprocess inside an mpi-multiprocess.

        Returns:
        :returns: bol bol:  Returns True if MPI can be run and there aren't any flags against running mpi mode
                        Returns False if it cannot or should not run mpi mode.
        """
        if MPI_installed == False:
            # mpi4py isn't installed and we will need to multiprocess
            return False

        if flag_for_low_level == True:
            # Flagged for low level and testing import mpi4py.MPI can be a problem
            return False

        if mode == "mpi" or mode == "None" or mode is None:
            # This must be either mpi or None, mpi4py can be installed and it hasn't been flagged at low level

            # Before executing Parallelizer with mpi4py (which override python raise Exceptions)
            # We must check that it is being run with the "-m mpi4py" runpy flag
            # Although this is lower priority over mpi4py version (as mpi4py.__versions__ less than 2.1.0 do not offer the -m feature)
            #      This should get checked before loading the mpi4py api
            if not mpi4py_launch_flag_present():
                print(MPI_LAUNCH_FLAG_MSG)
                return False

            try:
                return self.check_mpi_available()
            except Exception:
                return False

        return False

    # TODO Rename this here and in `test_import_MPI`
    def check_mpi_available(self) -> bool:
        """Report whether the installed mpi4py can actually be used for mpi.

        Importing the MPI sublibrary is the real test (it is what fails inside
        an already-mpi-parallelized program), so the import stays here even
        though nothing in this function uses the name.

        Returns:
            True if MPI imported and the mpi4py version is recent enough.
        """

        import mpi4py
        from mpi4py import MPI  # noqa: F401

        if not mpi4py_version_supported(mpi4py.__version__):
            print(MPI4PY_VERSION_MSG)
            return False

        return True

    def start(self, mode: str | None = None):
        """
        One must call this method before `run()` in order to configure MPI parallelization

        This creates the object for parallizing in a given mode.

        mode=None can be used at a top level program, but if using a program enabled
        with this multiprocess, referenced by a top level program, make sure the mode
        is explicitly chosen.

        Args:
            mode: the multiprocess mode to be used, ie) serial, multiprocessing, mpi, or None:
                if None then we will try to pick a possible multiprocessing choice. This should only be used for
                top level coding. It is best practice to specify which multiprocessing choice to use.
                if you have smaller programs used by a larger program, with both mpi enabled there will be problems, so specify multiprocessing is important.
        Returns:
        :returns: class parallel_obj: This is the obstantiated object of the class of parallizations.
                            ie) if self.mode=='mpi' self.parallel_obj will be an instance of the mpi class
                                This style of retained parallel_obj to be used later is important because this is the object which controls the work nodes and maintains the mpi universe
                            self.parallel_obj will be set to None for simpler parallization methods like serial
        """

        if mode is None:
            mode = self.mode

        if mode == "mpi":
            if self.HAS_MPI == True:
                # THIS IS EXPLICITILY CHOSEN TO BE RUN IN MPI AND CAN WORK WITH MPI
                ParallelMPI_obj = ParallelMPI()
                ParallelMPI_obj.start()
                return ParallelMPI_obj

            raise Exception("mpi4py package must be available to use mpi mode")

        return None

    def end(self, mode=None):
        """
        Call this method before exit to terminate MPI workers


        Inputs:
        :param str mode: the multiprocess mode to be used, ie) serial, multiprocessing, mpi, or None:
                            if None then we will try to pick a possible multiprocessing choice. This should only be used for
                            top level coding. It is best practice to specify which multiprocessing choice to use.
                            if you have smaller programs used by a larger program, with both mpi enabled there will be problems, so specify multiprocessing is important.
        """

        if mode is None:
            mode = self.mode
        if mode == "mpi":
            if self.HAS_MPI == True and self.parallel_obj != None:
                # THIS IS EXPLICITILY CHOSEN TO BE RUN IN MPI AND CAN WORK WITH MPI
                self.parallel_obj.end()

            else:
                raise Exception("mpi4py package must be available to use mpi mode")

    def run(self, args, func, num_procs=None, mode=None):
        """
        Run a task in parallel across the system.

        Mode can be one of 'mpi', 'multiprocessing' or 'none' (serial). If it is not
        set, the best value will be determined automatically.

        By default, this method will use the full resources of the system. However,
        if the mode is set to 'multiprocessing', num_procs can control the number
        of threads initialized when it is set to a nonzero value.

        Example: If one wants to multiprocess function  def foo(x,y) which takes 2 ints and one wants to test all permutations of x and y between 0 and 2:
                    args = [(0,0),(1,0),(2,0),(0,1),(1,1),(2,1),(0,2),(1,2),(2,2)]
                    func = foo      The namespace of foo


        Inputs:
        :param python_obj func: This is the object of the function which will be used.
        :param list args: a list of lists/tuples, each sublist/tuple must contain all information required by the function for a single object which will be multiprocessed
        :param int num_procs:  (Primarily for Developers)  the number of processors or nodes that will be used. If None than we will use all available nodes/processors
                                        This will be overriden and fixed to a single processor if mode==serial
        :param str mode:  (Primarily for Developers) the multiprocess mode to be used, ie) serial, multiprocessing, mpi, or None:
                            if None then we will try to pick a possible multiprocessing choice. This should only be used for
                            top level coding. It is best practice to specify which multiprocessing choice to use.
                            if you have smaller programs used by a larger program, with both mpi enabled there will be problems, so specify multiprocessing is important.
                            BEST TO LEAVE THIS BLANK
        Returns:
        :returns: list results: A list containing all the results from the multiprocess
        """

        # determine the mode
        if mode is None:
            mode = self.mode
        elif self.mode != mode:
            if mode not in ["mpi", "serial", "multiprocessing"]:
                printout = (
                    "Overriding function with a multiprocess mode which doesn't match: "
                    + mode
                )
                raise Exception(printout)
            if mode == "mpi":
                printout = "Overriding multiprocess can't go from non-mpi to mpi mode"
                raise Exception(printout)

        if num_procs is None:
            num_procs = self.num_procs

        if num_procs != self.num_procs and mode == "serial":
            printout = "Can't override num_procs in serial mode"
            raise Exception(printout)

        if mode != self.mode:
            printout = (
                f"changing mode from {self.mode} to {mode} for development purpose"
            )
            print(printout)

        # compute
        if mode == "mpi":
            if not self.HAS_MPI:
                raise Exception("mpi4py package must be available to use mpi mode")

            return self.parallel_obj.run(func, args)

        elif mode == "multiprocessing":
            return MultiThreading(args, num_procs, func)
        else:
            # serial is running the ParallelThreading with num_procs=1
            return MultiThreading(args, 1, func)

    def pick_mode(self):
        """
        Determines the parallelization cababilities of the system and returns one
        of the following modes depending on the configuration:

        Returns:
        :returns: str mode: the mode which is to be used 'mpi', 'multiprocessing', 'serial'
        """
        # check if mpi4py is loaded and we have more than one processor in our MPI world
        if self.HAS_MPI:
            try:
                if mpi4py.MPI.COMM_WORLD.Get_size() > 1:
                    return "mpi"
            except Exception:
                return "multiprocessing"
            else:
                return "multiprocessing"
        # # check if we could utilize more than one processor
        # if multiprocessing.cpu_count() > 1:
        #     return 'multiprocessing'

        # default to multiprocessing
        return "multiprocessing"

    def return_mode(self):
        """
        Returns the mode chosen for the parallelization cababilities of the system and returns one
        of the following modes depending on the configuration:
        :param str mode: the multiprocess mode to be used, ie) serial, multiprocessing, mpi, or None:
                    if None then we will try to pick a possible multiprocessing choice. This should only be used for
                    top level coding. It is best practice to specify which multiprocessing choice to use.
                    if you have smaller programs used by a larger program, with both mpi enabled there will be problems, so specify multiprocessing is important.
                    BEST TO LEAVE THIS BLANK
        Returns:
        :returns: str mode: the mode which is to be used 'mpi', 'multiprocessing', 'serial'
        """
        return self.mode

    def compute_nodes(self, mode=None):
        """
        Computes the number of "compute nodes" according to the selected mode.

        For mpi, this is the universe size
        For multiprocessing this is the number of available cores
        For serial, this value is 1
        Returns:
        :returns: int num_procs: the number of nodes/processors which is to be used
        """
        if mode is None:
            mode = self.mode

        if mode == "mpi":
            if not self.HAS_MPI:
                raise Exception("mpi4py package must be available to use mpi mode")
            return mpi4py.MPI.COMM_WORLD.Get_size()
        elif mode == "multiprocessing":
            return multiprocessing.cpu_count()
        else:
            return 1

    def return_node(self):
        """
        Returns the number of "compute nodes" according to the selected mode.

        For mpi, this is the universe size
        For multiprocessing this is the number of available cores
        For serial, this value is 1
        Returns:
        :returns: int num_procs: the number of nodes/processors which is to be used
        """
        return self.num_procs


class ParallelMPI(object):
    """
    Utility code for running tasks in parallel across an MPI cluster.
    """

    def __init__(self):
        """
        Default num_procs is all the processesors possible
        """

        self.COMM = mpi4py.MPI.COMM_WORLD

        self.Empty_object = Empty_obj()

    def start(self):
        """
        Call this method at the beginning of program execution to put non-root processors
        into worker mode.
        """

        rank = self.COMM.Get_rank()

        if rank == 0:
            return
        else:
            worker = self._worker()

    def end(self):
        """
        Call this method to terminate worker processes
        """

        self.COMM.bcast(None, root=0)

    def _worker(self):
        """
        Worker processors wait in this function to receive new jobs
        """
        while True:
            # receive function for new job
            func = self.COMM.bcast(None, root=0)

            # kill signal
            if func is None:
                exit(0)

            # receive arguments, then the seeds the root drew for them. Both
            # scatters run on every rank in the same order, whether or not
            # this rank got real work; otherwise the collectives fall out of
            # step.
            args_chunk = self.COMM.scatter([], root=0)
            seed_chunk = self.COMM.scatter([], root=0)

            if type(args_chunk[0]) == type(
                self.Empty_object
            ):  # or  args_chunk[0] == [[self.Empty_object]]:
                result_chunk = [[self.Empty_object]]
                result_chunk = self.COMM.gather(result_chunk, root=0)

            else:
                # perform the calculation and send results. run_one turns a
                # raised exception into None, matching every other dispatch
                # path; calling func bare here let one bad job skip the gather
                # below, which blocks every other rank forever.
                result_chunk = [
                    run_one(func, arg, seed[0])
                    for arg, seed in zip(args_chunk, seed_chunk)
                    if type(arg[0]) != type(self.Empty_object)
                ]
                result_chunk = self.COMM.gather(result_chunk, root=0)

    def handle_undersized_jobs(self, arr, n):
        if len(arr) > n:
            printout = "the length of the package is bigger than the length of the number of nodes!"
            print(printout)
            raise Exception(printout)

        filler_slot = [[self.Empty_object]]
        while len(arr) < n:
            arr.append(filler_slot)
            if len(arr) == n:
                break

        return arr

    def _split(self, arr, n):
        """
        Takes an array of items and splits it into n "equally" sized
        chunks that can be provided to a worker cluster.
        """

        s = len(arr) // n
        remainder = len(arr) - int(s) * int(n)

        chuck_list = []
        temp = []
        counter = 0
        for x in range(len(arr)):
            # add 1 per group until remainder is removed
            if remainder != 0:
                r = 1
                if counter == s + 1:
                    remainder = remainder - 1
            else:
                r = 0

            if counter == s + r:
                chuck_list.append(temp)
                temp = []
                counter = 1
            else:
                counter += 1

            temp.append(list(arr[x]))
            if x == len(arr) - 1:
                chuck_list.append(temp)

        if len(chuck_list) != n:
            chuck_list = self.handle_undersized_jobs(chuck_list, n)

        return chuck_list

    def _join(self, arr):
        """
        Joins a "list of lists" that was previously split by _split().

        Returns a single list.
        """
        arr = tuple(arr)
        arr = [x for x in arr if type(x) != type(self.Empty_object)]
        arr = [a for sub in arr for a in sub]
        return [x for x in arr if type(x) != type(self.Empty_object)]

    def check_and_format_args(self, args):
        # Make sure args is a list of lists
        if type(args) not in [list, tuple]:
            printout = "args must be a list of lists"
            print(printout)
            raise Exception(printout)

        item_type = type(args[0])
        for i in range(len(args)):
            if type(args[i]) == item_type:
                continue
            printout = "all items within args must be the same type and must be either a list or tuple"
            print(printout)
            raise Exception(printout)
        if item_type == list:
            return args
        elif item_type == tuple:
            args = [list(x) for x in args]
            return args
        else:
            printout = "all items within args must be either a list or tuple"
            print(printout)
            raise Exception(printout)

    def run(self, func, args):
        """
        Run a function in parallel across the current MPI cluster.

        * func is a pure function of type (A)->(B)
        * args is a list of type list(A)

        This method batches the computation across the MPI cluster and returns
        the result of type list(B) where result[i] = func(args[i]).

        Important note: func must exist in the namespace at initialization.
        """
        num_of_args_start = len(args)
        if len(args) == 0:
            return []
        args = self.check_and_format_args(args)

        size = self.COMM.Get_size()

        # One seed per job, drawn here in job order. _split depends only on the
        # list length and the rank count, so chunking the seeds the same way
        # lands each job's seed on the rank that got the job.
        seeds = [[seed] for seed in draw_job_seeds(len(args))]

        # broadcast function to worker processors
        self.COMM.bcast(func, root=0)

        # chunkify the argument list
        args_chunk = self._split(args, size)
        seed_chunks = self._split(seeds, size)

        # scatter argument chunks to workers
        args_chunk = self.COMM.scatter(args_chunk, root=0)
        seed_chunk = self.COMM.scatter(seed_chunks, root=0)

        if type(args_chunk) != list:
            raise Exception("args_chunk needs to be a list")

        # perform the calculation and get results
        result_chunk = [
            run_one(func, arg, seed[0]) for arg, seed in zip(args_chunk, seed_chunk)
        ]
        sys.stdout.flush()

        result_chunk = self.COMM.gather(result_chunk, root=0)

        if type(result_chunk) != list:
            raise Exception("result_chunk needs to be a list")

        # group results
        results = self._join(result_chunk)

        if len(results) != num_of_args_start:
            results = [x for x in results if type(x) != type(self.Empty_object)]
            results = flatten_list(results)

        if type(results) != list:
            raise Exception("results needs to be a list")

        results = [x for x in results if type(x) != type(self.Empty_object)]
        sys.stdout.flush()
        return results


#


class Empty_obj(object):
    """
    Create a unique Empty Object to hand to empty processors
    """

    pass


#


"""
Run commands on multiple processors in python.

Adapted from examples on https://docs.python.org/2/library/multiprocessing.html
"""


def seed_generators(seed: int) -> None:
    """Seed both global generators Gypsum-DL samples from, from one value.

    Variant selection samples with the random module and the ring-conformation
    clustering samples with numpy (by way of scipy's kmeans2), so seeding one
    without the other still leaves the caller drawing from an unseeded stream.
    Both are derived from a single value so that every seeded context in the
    program (the parent process at startup, each individual job) sets up its
    streams the same way.

    Args:
        seed: The seed. numpy rejects seeds that do not fit in 32 bits, while
            the random module accepts an integer of any size, so larger values
            are folded in rather than failing the run over the choice of seed.
    """

    random.seed(seed)
    numpy.random.seed(seed % 2**32)


def draw_job_seeds(num_jobs: int) -> list[int]:
    """Draw one seed per job from the caller's generator, in job order.

    Seeding per worker rather than per job cannot make a run reproducible: the
    task queue is drained by whichever worker is free, so a job's stream would
    depend on who picked it up. Drawing here instead, in the dispatching
    process before any job runs, ties a job's stream to its position in the
    input alone, which makes a seeded run reproducible whatever the job manager
    and worker count. An unseeded parent still yields distinct per-job seeds,
    which is what keeps forked workers off a shared numpy stream.

    Args:
        num_jobs: How many seeds to draw.

    Returns:
        One 32-bit seed per job, positionally matching the job list.
    """

    return [random.getrandbits(32) for _ in range(num_jobs)]


def run_one(
    func: Callable[..., _JobResult],
    args: Sequence[object],
    seed: int | None = None,
) -> _JobResult | None:
    """Run a single job, turning a raised exception into a None result.

    Every non-MPI dispatch path routes through here (the multiprocessing
    worker, the single-processor branch of MultiThreading, and the in-process
    branch of each pipeline step) so that a molecule that raises is dropped
    identically no matter which job_manager is in use. Callers already treat
    None as a failed job and strip it.

    Args:
        func: The job function to call.
        args: Positional arguments to unpack into `func`.
        seed: Seed for this job's generators, or None to draw from the state
            already in place. The caller's own generator state is restored
            afterwards, so a job run in-process (serial mode, the mpi root's
            chunk) leaves the dispatching stream exactly where a job run in a
            child process would.

    Returns:
        Whatever `func` returns, or None if `func` raised.
    """
    if seed is not None:
        saved_random_state = random.getstate()
        saved_numpy_state = numpy.random.get_state()
        seed_generators(seed)

    try:
        return func(*args)
    except Exception:
        name = getattr(func, "__name__", repr(func))
        print(f"ERROR in {name}: {traceback.format_exc()}")
        return None
    finally:
        if seed is not None:
            random.setstate(saved_random_state)
            numpy.random.set_state(saved_numpy_state)


def MultiThreading(inputs, num_procs, task_name):
    """Initialize this object.

    Args:
        inputs ([data]): A list of data. Each datum contains the details to
            run a single job on a single processor.
        num_procs (int): The number of processors to use.
        task_class_name (class): The class that governs what to do for each
            job on each processor.
    """

    results = []

    # If there are no inputs, just return an empty list.
    if len(inputs) == 0:
        return results

    inputs = check_and_format_inputs_to_list_of_tuples(inputs)

    num_procs = count_processors(len(inputs), num_procs)

    tasks = []

    # Every seed is drawn here, before any job runs, so that job k gets the
    # same stream whether it runs in this process or in whichever worker
    # happens to pull it off the queue.
    seeds = draw_job_seeds(len(inputs))

    for index, item in enumerate(inputs):
        if not isinstance(item, tuple):
            item = (item,)
        task = (index, (task_name, item, seeds[index]))
        tasks.append(task)

    if num_procs == 1:
        for item in tasks:
            job, args, seed = item[1]
            results.append(run_one(job, args, seed))
    else:
        results = start_processes(tasks, num_procs)

    return results


###
# Worker function
###


def worker(
    input: "multiprocessing.Queue[object]",
    output: "multiprocessing.Queue[object]",
) -> None:
    """Consume jobs from a queue and report their results.

    Each job carries the seed drawn for it by the dispatching process, and
    run_one installs that seed before calling the job. That matters twice over.
    Under the fork start method a child inherits a byte-identical copy of the
    parent's numpy legacy global state and nothing reseeds it (unlike the
    random module, which CPython reseeds in the child through
    os.register_at_fork, numpy installs no such hook), so without a reseed
    every worker would cluster ring conformations against the same draw. And
    because the seed travels with the job rather than with the worker, the
    result does not depend on which worker drained which task.

    Args:
        input: Queue of (index, (function, arguments, seed)) jobs, terminated
            by the string "STOP".
        output: Queue the (index, result) pairs are reported on.
    """

    for seq, job in iter(input.get, "STOP"):
        func, args, seed = job
        # A dead worker would leave the parent blocked on done_queue.get()
        # forever, so failures must be reported as results.
        output.put((seq, run_one(func, args, seed)))


def check_and_format_inputs_to_list_of_tuples(args):
    # Make sure args is a list of tuples
    if type(args) not in [list, tuple]:
        printout = "args must be a list of tuples"
        print(printout)
        raise Exception(printout)

    item_type = type(args[0])
    for i in range(len(args)):
        if type(args[i]) == item_type:
            continue

        printout = "all items within args must be the same type and must be either a list or tuple"
        print(printout)
        raise Exception(printout)
    if item_type == tuple:
        return args
    elif item_type == list:
        args = [tuple(x) for x in args]
        return args
    else:
        printout = "all items within args must be either a list or tuple"
        print(printout)
        raise Exception(printout)


def count_processors(num_inputs, num_procs):
    """
    Checks processors available and returns a safe number of them to
    utilize.

    :param int num_inputs: The number of inputs.
    :param int num_procs: The number of desired processors.

    :returns: The number of processors to use.
    """
    # first, if num_procs <= 0, determine the number of processors to
    # use programatically
    if num_procs <= 0:
        num_procs = multiprocessing.cpu_count()

    # reduce the number of processors if too many have been specified
    if num_inputs < num_procs:
        num_procs = num_inputs

    return num_procs


def start_processes(
    inputs: Sequence[
        tuple[int, tuple[Callable[..., _JobResult], Sequence[object], int]]
    ],
    num_procs: int,
) -> list[_JobResult | None]:
    """Run the given jobs across worker processes and collect their results.

    Results come back in input order. A worker that dies outright (the OOM
    killer, a native fault in RDKit) never reports the result it was holding,
    so the collection loop is bounded rather than blocking forever on a result
    that will never arrive: it waits in intervals and gives up once no worker
    is left to produce anything. The workers are also joined before returning,
    which keeps finished runs from leaving non-daemon children behind.

    Args:
        inputs: Jobs as (index, (function, arguments, seed)) pairs, where the
            index determines the position of the job's result in the returned
            list and the seed is installed by the worker before the job runs.
        num_procs: How many worker processes to start.

    Returns:
        One entry per input, in input order; None where the job raised.

    Raises:
        Exception: If every worker exited before all results were reported.
    """

    # Create queues
    task_queue: multiprocessing.Queue[object] = multiprocessing.Queue()
    done_queue: multiprocessing.Queue[tuple[int, _JobResult | None]] = (
        multiprocessing.Queue()
    )

    # Submit tasks
    for item in inputs:
        task_queue.put(item)

    # Queue the stop sentinels up front. The queue is FIFO and every real task
    # is already enqueued, so no worker can see a sentinel early; doing it here
    # rather than after the collection loop means the workers still shut down
    # if the loop exits by raising.
    for _ in range(num_procs):
        task_queue.put("STOP")

    # Start worker processes. Children inherit the parent's numpy generator
    # state on fork; the per-job seed each task carries is what keeps them off
    # that shared stream.
    procs = [
        multiprocessing.Process(target=worker, args=(task_queue, done_queue))
        for _ in range(num_procs)
    ]
    for proc in procs:
        proc.start()

    results: list[tuple[int, _JobResult | None]] = []
    while len(results) < len(inputs):
        try:
            results.append(done_queue.get(timeout=WORKER_POLL_TIMEOUT))
            continue
        except queue.Empty:
            pass

        if any(proc.is_alive() for proc in procs):
            continue

        # No worker is running. A result put just before the last worker exited
        # can still be in transit, so drain once more before concluding that
        # one was lost.
        try:
            results.append(done_queue.get(timeout=WORKER_EXIT_GRACE))
            continue
        except queue.Empty:
            pass

        for proc in procs:
            proc.join()

        # Unconsumed tasks may still be buffered for the feeder thread, which
        # would block interpreter exit now that nothing is reading the queue.
        task_queue.cancel_join_thread()

        raise Exception(
            f"A worker process died after returning {len(results)} of "
            f"{len(inputs)} results."
        )

    for proc in procs:
        proc.join()

    results.sort(key=lambda tup: tup[0])

    return [item[1] for item in map(list, results)]


def flatten_list(tier_list: list) -> list:
    """
    Given a list of lists, this returns a flat list of all items.

    Args:
        tier_list: A 2D list.

    Returns:
        A flat list of all items.
    """

    if tier_list is None:
        return []
    flattened: list = []
    for item in tier_list:
        if item is None:
            continue
        if isinstance(item, list):
            flattened.extend(x for x in item if x is not None)
        else:
            flattened.append(item)
    return flattened


def strip_none(none_list: list[Any]) -> list[Any]:
    """
    Given a list that might contain None items, this returns a list with no
    None items.

    :params list none_list: A list that may contain None items.

    :returns: A list stripped of None items.
    """

    return [] if none_list is None else [x for x in none_list if x is not None]
