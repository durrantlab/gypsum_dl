import sys
import argparse
from gypsum.Start import ConfGenerator

parser = argparse.ArgumentParser()
parser.add_argument('--json', '-j',
                    help='name of a json file containing all parameters')
parser.add_argument('--input', '-i',
                    help='name of the input file')
parser.add_argument('--output', '-o',
                    help='name of the output file')
parser.add_argument('--source', '-s',
                    help='TO BE DETERMINED')
parser.add_argument('--separate_output_files', action='store_true',
                    help='indicates that the output should be split')
parser.add_argument('--openbabel_exe',
                    help='path to the openbabel executable')
parser.add_argument('--processors', '-p', help='number of processors to use in parallel')
parser.add_argument('--min_ph', help='minimum pH to consider')
parser.add_argument('--max_ph', help='maximum pH to consider')
parser.add_argument('--delta_ph', help='')
parser.add_argument('--thoroughness', '-t',
                    help='')
parser.add_argument('--max_variants', '-m',
                    help='')
parser.add_argument('--optimize_geometry', action='store_true',
                    help='')
parser.add_argument('--alternate_ring_conformations', action='store_true',
                    help='')

args = parser.parse_args()
# # Test chiral. One chiral center specified, one not.
# ConfGenerator(['I[C@@](C)(F)C(I)(F)C'], 5, 9)

# # Test double-bond enumeration. Should be two.
# ConfGenerator(['FC(I)=C(Cl)C'], 5, 9)

# # Test Tautomers
# ConfGenerator(['C1=CC=CNC(=O)1'], 5, 9)

# # Test Tautomers from including bad tautomers. This also tests protonation at various states.
# ConfGenerator(['N[C@@H](CC1=CNC=N1)C(O)=O'], 5, 9)

# # Test non-aromatic rings. Should generate multiple ring confomers.
#ConfGenerator(['C1(C=CC=C2)=C2CCCC1'], 5, 9)

# Test all
#ConfGenerator(['C1(C=CC=C2)=C2CCCC1', 'N[C@@H](CC1=CNC=N1)C(O)=O', 'C1=CC=CNC(=O)1', 'FC(I)=C(Cl)C', 'I[C@@](C)(F)C(I)(F)C'], 5, 9)
#ConfGenerator("samples.smi", 5, 9)

#ConfGenerator("sample_parameters.json")
ConfGenerator(sys.argv[1])
