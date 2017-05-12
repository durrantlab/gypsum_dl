#!/usr/bin/env python
"""
Gypsum is a conversion script to transform smiles strings and 2D SDFs
into 3D models.
"""
import argparse
import copy
from gypsum.Start import ConfGenerator


PARSER = argparse.ArgumentParser()


PARSER.add_argument('--json', '-j', metavar='param.json',
                    help='Name of a json file containing all parameters. \
                    Overrides other arguments.')
PARSER.add_argument('--source', '-s', metavar='input.smi',
                    help='Name of the source file.')
PARSER.add_argument('--output_file', '-o', metavar='output.sdf',
                    help='Name of the output file.')
PARSER.add_argument('--openbabel_executable',
                    help='Path to the openbabel executable.')
PARSER.add_argument('--num_processors', '-p', metavar='N', default=1,
                    help='Number of processors to use in parallel.')
PARSER.add_argument('--min_ph', metavar='MIN', type=float,
                    help='Minimum pH to consider.')
PARSER.add_argument('--max_ph', metavar='MAX', type=float,
                    help='Maximum pH to consider.')
PARSER.add_argument('--delta_ph_increment', metavar='D', type=float,
                    help='Increments of pH between min and max.')
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

ARGS_DICT = vars(PARSER.parse_args())

INPUTS = copy.deepcopy(ARGS_DICT)

for k, v in ARGS_DICT.items():
    if v is None:
        del INPUTS[k]

GENERATOR = ConfGenerator(INPUTS)
GENERATOR.run()
