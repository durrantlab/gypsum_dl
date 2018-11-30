#!/usr/bin/env python
"""
Gypsum is a conversion script to transform smiles strings and 2D SDFs
into 3D models.
"""
import argparse
import copy
def handle_mpi(ARGS_DICT):

    MPI_installed = False
    try:
        import mpi4py 
        import mpi4py.MPI 
        MPI_installed = True
    except:
        MPI_installed = False

    if MPI_installed == False:
        printout = "Tried to be run in MPI mode without mpi4py dependency. Please install mpi4py or remove the --mpi from user params"
        print(printout)
        raise Exception(printout)
    
    COMM = mpi4py.MPI.COMM_WORLD

    rank = COMM.Get_rank()

    if rank == 0:
        params = get_params(ARGS_DICT)
    else:
        pass
    COMM.barrier()
    from gypsum.Start import conf_generator
    conf_generator(params)  
    print("Finished Gypsum")


def get_params(ARGS_DICT):
    INPUTS = copy.deepcopy(ARGS_DICT)

    for k, v in ARGS_DICT.items():
        if v is None:
            del INPUTS[k]

    from gypsum.Start import get_user_params
    params = get_user_params(INPUTS)

    return params

def run_gypsum(ARGS_DICT):

    params = get_params(ARGS_DICT)
    from gypsum.Start import conf_generator
    conf_generator(params)  
    print("Finished Gypsum")

#

PARSER = argparse.ArgumentParser()


PARSER.add_argument('--mpi', '-q', action='store_true',
                    help='place in bash command ex) python run_gypsum --mpi -j json_file.json')
PARSER.add_argument('--json', '-j', metavar='param.json',
                    help='Name of a json file containing all parameters. \
                    Overrides other arguments.')
PARSER.add_argument('--source', '-s', metavar='input.smi',
                    help='Name of the source file.')
PARSER.add_argument('--output_file', '-o', metavar='output.sdf',
                    help='Name of the output file.')
PARSER.add_argument('--num_processors', '-p', metavar='N', default=1,
                    help='Number of processors to use in parallel.')
PARSER.add_argument('--min_ph', metavar='MIN', type=float,
                    help='Minimum pH to consider.')
PARSER.add_argument('--max_ph', metavar='MAX', type=float,
                    help='Maximum pH to consider.')
PARSER.add_argument('--ph_std_dev', metavar='D', type=float,
                    help='Size of pH substructure ranges.')
PARSER.add_argument('--thoroughness', '-t',
                    help='How wide a search to look for conformers.')
PARSER.add_argument('--max_variants_per_compound', '-m', type=int, metavar='V',
                    help='The maximum number of models to create.')
PARSER.add_argument('--separate_output_files', action='store_true',
                    help='Indicates that the outputs should be split between files.')
PARSER.add_argument('--skip_optimize_geometry', action='store_true',
                    help='Skips the optimization step.')
PARSER.add_argument('--skip_alternate_ring_conformations', action='store_true',
                    help='Skips using non-aromatic ring conformations.')
PARSER.add_argument('--skip_adding_hydrogen', action='store_true',
                    help='Skips adding hydrogens based on pH.')
PARSER.add_argument('--skip_making_tautomers', action='store_true',
                    help='Skips tautomer creation.')
PARSER.add_argument('--skip_ennumerate_chiral_mol', action='store_true',
                    help='Skips the ennumeration of chiral centers.')
PARSER.add_argument('--skip_ennumerate_double_bonds', action='store_true',
                    help='Skips the ennumeration of double bonds.')
PARSER.add_argument('--2d_output_only', action='store_true',
                    help='Skips using non-aromatic ring conformations.')
PARSER.add_argument('--multithread_mode', default='multithreading', choices = ["mpi","multithreading","serial"],
                    help='Determine what style multithreading: mpi, multithreading, or serial.\
                    If this program is being used by a program in MPI mode we recommend setting this to serial.\
                    serial will override num_processors and force it to be on a single processor.')

ARGS_DICT = vars(PARSER.parse_args())

if ARGS_DICT["mpi"] == True:
    handle_mpi(ARGS_DICT)
else:
    run_gypsum(ARGS_DICT)
import sys
sys.exit(0)



