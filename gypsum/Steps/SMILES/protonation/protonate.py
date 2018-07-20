"""
Main entry script for running protonation
"""
import argparse

import protonation_functions as pf

parser = argparse.ArgumentParser(description="Protonates small moleucles.")
parser.add_argument('--min_ph', metavar='MIN', type=float, default=6.4,
                    help='Minimum pH to consider.')
parser.add_argument('--max_ph', metavar='MAX', type=float, default=8.4,
                    help='Maximum pH to consider.')
parser.add_argument('--st_dev', metavar='STD', type=float, default=1.0,
                    help='Standard devation range (number of standard devs).')
parser.add_argument('--smiles', type=str,
                    help='SMILE string to protonate.')
parser.add_argument('--smiles_file', type=str,
                    help='File which contains SMILES strings to protonate.')
parser.add_argument('--output_file', type=str,
                    help='File to write protonated SMILES. (Optional)')


if __name__ == "__main__":
    args = vars(parser.parse_args())

    output = pf.protonate(args)

    if 'output_file' in args:
        with open(args['output_file'], 'w') as file:
            for out in output:
                file.write("\t".join(out))
    else:
        for out in output:
            print(out)
