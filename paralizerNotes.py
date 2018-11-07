# # 
# old style
import gypsum.gypsum.Multiprocess as mp

job_input = [(mgl_python, receptor_template, i)  for i in list_of_jobs]

# Fix this to 1 processor, so no overwriting issues, but could be expanded if we ever wanted to do multiple receptors
output = mp.MultiThreading(job_input, 1,  foo)
    #


####################################################
# New usage


#At start code one must create the parallelizer object which will be passaged through the full script

from gypsum.parallelizer import Parallelizer


# # # launch mpi workers
if params["multithread_mode"] == 'mpi':
    parallelizer_obj = Parallelizer(params["multithread_mode"], params["num_processors"])
else:
    # Lower level mpi (ie making a new Parallelizer within an mpi) 
    #   has problems with importing the MPI enviorment and mpi4py
    #   So we will flag it to skip the MPI mode and just go to multithread/serial
    # This is a saftey precaution
    parallelizer_obj = Parallelizer(params["multithread_mode"], params["num_processors"], True)


# parallelizer_obj MUST BE THREADED THROUGHOUT OR MPI WILL BREAK BADLY

# parallelizer_obj contains information about the processor type and number of node/processors in self
#   So its unnecessary to provide multithread_mode or num_processors again but it is possible


    tmp = Parallelizer_obj.run(parallel_addH, inputs, num_processors, multithread_mode)

    tmp = Parallelizer_obj.run(parallel_addH, inputs)