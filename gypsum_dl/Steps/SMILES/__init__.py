import os

# gypsum_dl/gypsum_dl/Steps/SMILES/
# Including the below allows other programs to import functions from
# gypsum-DL.
import sys

current_dir = os.path.dirname(os.path.realpath(__file__))
Steps = os.path.dirname(current_dir)
gypsum_gypsum_dir = os.path.dirname(Steps)
gypsum_top_dir = os.path.dirname(gypsum_gypsum_dir)
sys.path.extend([current_dir, Steps, gypsum_gypsum_dir, gypsum_top_dir])
