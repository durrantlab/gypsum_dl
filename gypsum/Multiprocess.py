"""
Run commands on multiple processors in python.

Adapted from examples on https://docs.python.org/2/library/multiprocessing.html
"""
import __future__


import multiprocessing


def MultiThreading(inputs, num_processors, task_name):
    """Initialize this object.

    Args:
        inputs ([data]): A list of data. Each datum contains the details to
            run a single job on a single processor.
        num_processors (int): The number of processors to use.
        task_class_name (class): The class that governs what to do for each
            job on each processor.
    """

    results = []

    # If there are no inputs, just return an empty list.
    if len(inputs) == 0:
        return results

    num_processors = count_processors(len(inputs), num_processors)

    tasks = []

    for index, item in enumerate(inputs):
        if not isinstance(item, tuple):
            item = (item,)
        task = (index, (task_name, item))
        tasks.append(task)

    if num_processors == 1:
        for item in tasks:
            job, args = item[1]
            output = job(*args)
            results.append(output)
    else:
        results = start_processes(tasks, num_processors)


###
# Worker function
###

def worker(input, output):
    for seq, job in iter(input.get, 'STOP'):
        func, args = job
        result = func(*args)
        ret_val = (seq, result)
        output.put(ret_val)



def count_processors(num_inputs, num_processors):
    """
    Checks processors available and returns a safe number of them to
    utilize.

    :param int num_inputs: The number of inputs.
    :param int num_processors: The number of desired processors.

    :returns: The number of processors to use.
    """
    # first, if num_processors <= 0, determine the number of processors to
    # use programatically
    if num_processors <= 0:
        num_processors = multiprocessing.cpu_count()

    # reduce the number of processors if too many have been specified
    if num_inputs < num_processors:
        num_processors = num_inputs

    return num_processors

def start_processes(inputs, num_processors):
    """
    Creates a queue of inputs and outputs
    """

    # Create queues
    task_queue = multiprocessing.Queue()
    done_queue = multiprocessing.Queue()

    # Submit tasks
    for item in inputs:
        task_queue.put(item)

    # Start worker processes
    for i in range(num_processors):
        multiprocessing.Process(target=worker, args=(task_queue, done_queue)).start()

    # Get and print results
    results = []
    for i in range(len(inputs)):
        results.append(done_queue.get())

    # Tell child processes to stop
    for i in range(num_processors):
        task_queue.put('STOP')

    results.sort(key=lambda tup: tup[0])

    return  [item[1] for item in map(list, results)]

###
# Helper functions
###

def flatten_list(tier_list):
    """
    Given a list of lists, this returns a flat list of all items.

    :params list tier_list: A 2D list.

    :returns: A flat list of all items.
    """
    if tier_list is None:
        return []
    flat_list = [item for sublist in tier_list for item in sublist]
    return flat_list

def strip_none(none_list):
    """
    Given a list that might contain None items, this returns a list with no
    None items.

    :params list none_list: A list that may contain None items.

    :returns: A list stripped of None items.
    """
    if none_list is None:
        return []
    results = [x for x in none_list if x != None]
    return results

