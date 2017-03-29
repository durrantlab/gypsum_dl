"""
Run commands on multiple processors in python.

Adapted from examples on https://docs.python.org/2/library/multiprocessing.html
"""



import multiprocessing



def MultiThreading(inputs, num_processors, task_class_name):
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

    for item in inputs:
        if not isinstance(item, tuple):
            item = (item,)
        tasks.append((task_class_name, item))

    results = start_processes(tasks, num_processors)

    return results

#
# Worker function
#
def worker(input, output):
    for func, args in iter(input.get, 'STOP'):
        result = func(*args)
        output.put(result)

#
# 
#
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

    return results
