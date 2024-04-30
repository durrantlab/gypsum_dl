import os

# gypsum_dl/gypsum_dl/
# Including the below allows other programs to import functions from
# gypsum-DL.
import sys

current_dir = os.path.dirname(os.path.realpath(__file__))
gypsum_gypsum_dir = current_dir
gypsum_top_dir = os.path.dirname(gypsum_gypsum_dir)
sys.path.extend([gypsum_gypsum_dir, gypsum_top_dir])
