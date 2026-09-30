"""Run Gypsum-DL as a module.

The installed gypsum-dl console script cannot be handed to mpi4py by module
name, so this is what makes "python -m mpi4py -m gypsum_dl" work under mpirun.
"""

from gypsum_dl.run import main

if __name__ == "__main__":
    main()
