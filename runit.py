import sys
import argparse
from gypsum.Start import ConfGenerator

parser = argparse.ArgumentParser()
parser.add_argument('--json', '-j', metavar='param.json',
                    help='Name of a json file containing all parameters. Overrides other arguments.')
parser.add_argument('--source', '-s', metavar='input.smi',
                    help='Name of the source file.')
parser.add_argument('--output_file', '-o', metavar='output.sdf',
                    help='Name of the output file.')
parser.add_argument('--openbabel_executable',
                    help='Path to the openbabel executable.')
parser.add_argument('--num_processors', '-p', metavar='N', default=1,
                    help='Number of processors to use in parallel.')
parser.add_argument('--min_ph', metavar='MIN', type=float,
                    help='Minimum pH to consider.')
parser.add_argument('--max_ph', metavar='MAX', type=float,
                    help='Maximum pH to consider.')
parser.add_argument('--delta_ph_increment', metavar='D', type=float,
                    help='Increments of pH between min and max.')
parser.add_argument('--thoroughness', '-t',
                    help='How wide a search to look for conformers.')
parser.add_argument('--max_variants_per_compound', '-m', type=int, metavar='V',
                    help='The maximum number of models to create.')
parser.add_argument('--separate_output_files', action='store_true',
                    help='Indicates that the outputs should be split between files.')

parser.add_argument('--skip_optimize_geometry', action='store_true',
                    help='Skips the optimization step.')
parser.add_argument('--skip_alternate_ring_conformations', action='store_true',
                    help='Skips using non-aromatic ring conformations.')

parser.add_argument('--skip_adding_hydrogen', action='store_true',
                    help='Skips adding hydrogens based on pH.')
parser.add_argument('--skip_making_tautomers', action='store_true',
                    help='Skips tautomer creation.')
parser.add_argument('--skip_ennumerate_chiral_mol', action='store_true',
                    help='Skips .')
parser.add_argument('--skip_ennumerate_double_bonds', action='store_true',
                    help='Skips .')

parser.add_argument('--2d_output_only', action='store_true',
                    help='Skips using non-aromatic ring conformations.')

args_dict = vars(parser.parse_args())

for k, v in args_dict.items():
    if v is None:
        del args_dict[k]

ConfGenerator(args_dict)
