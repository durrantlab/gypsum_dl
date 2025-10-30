import argparse
import copy

from gypsum_dl import utils
from gypsum_dl.start import prepare_molecules


def main():
    PARSER = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="""
    Gypsum-DL 1.2.2, a free, open-source program for preparing 3D small-molecule
    models. Beyond simply assigning atomic coordinates, Gypsum-DL accounts for
    alternate ionization, tautomeric, chiral, cis/trans isomeric, and
    ring-conformational forms.""",
        epilog="""
    EXAMPLES OF USE:

    1. Prepare a virtual library and save all 3D models to a single SDF file in the
    present directory:

    python run_gypsum_dl.py --source ./examples/sample_molecules.smi

    2. Instead save all 3D models to a different, existing folder:

    python run_gypsum_dl.py --source ./examples/sample_molecules.smi \\
    --output_folder /my/folder/

    3. Additionally save the models associated with each input molecule to
    separate files:

    python run_gypsum_dl.py --source ./examples/sample_molecules.smi \\
        --output_folder /my/folder/ --separate_output_files

    4. In addition to saving a 3D SDF file, also save 3D PDB files and an HTML file
    with 2D structures (for debugging).

    python run_gypsum_dl.py --source ./examples/sample_molecules.smi \\
        --output_folder /my/folder/ --add_pdb_output --add_html_output

    5. Save at most two variants per input molecule:

    python run_gypsum_dl.py --source ./examples/sample_molecules.smi \\
        --output_folder /my/folder/ --max_variants_per_compound 2

    6. Control how Gypsum-DL ionizes the input molecules:

    python run_gypsum_dl.py --source ./examples/sample_molecules.smi \\
        --output_folder /my/folder/ --min_ph 12 --max_ph 14 --pka_precision 1

    7. Run Gypsum-DL in serial mode (using only one processor):

    python run_gypsum_dl.py --source ./examples/sample_molecules.smi \\
        --job_manager serial

    8. Run Gypsum-DL in multiprocessing mode, using 4 processors:

    python run_gypsum_dl.py --source ./examples/sample_molecules.smi \\
        --job_manager multiprocessing --num_processors 4

    9. Run Gypsum-DL in mpi mode using all available processors:

    mpirun -n $NTASKS python -m mpi4py run_gypsum_dl.py \\
        --source ./examples/sample_molecules.smi \\
        --job_manager mpi --num_processors -1

    10. Gypsum-DL can also take parameters from a JSON file:

    python run_gypsum_dl.py --json myparams.json

    Where myparams.json might look like:

    {
        "source": "./examples/sample_molecules.smi",
        "separate_output_files": true,
        "job_manager": "multiprocessing",
        "output_folder": "/my/folder/",
        "add_pdb_output": true,
        "add_html_output": true,
        "num_processors": -1
    }
    """,
    )

    PARSER.add_argument(
        "--json",
        "-j",
        type=str,
        metavar="param.json",
        help="Name of a json file containing all parameters. \
                        Overrides all other arguments specified at the commandline.",
    )
    PARSER.add_argument(
        "--source",
        "-s",
        type=str,
        metavar="input.smi",
        help="Name of the source file (e.g., input.smi). Note: support for SMI (SMILES) files is better than support for SDF files, though Gypsum-DL can handle both.",
    )
    PARSER.add_argument(
        "--output_folder",
        "-o",
        type=str,
        help="The path to an existing folder where the Gypsum-DL "
        + "output file(s) will be saved.",
    )
    PARSER.add_argument(
        "--job_manager",
        type=str,
        default="multiprocessing",
        choices=["mpi", "multiprocessing", "serial"],
        help="Determine what style of multiprocessing to use: mpi, \
                            multiprocessing, or serial. Serial will override the \
                            num_processors flag, forcing it to be one. MPI mode \
                            requires mpi4py 2.1.0 or higher and should be executed \
                            as: mpirun -n $NTASKS python -m mpi4py run_gypsum_dl.py \
                            ...-settings...",
    )
    PARSER.add_argument(
        "--num_processors",
        "-p",
        type=int,
        metavar="N",
        default=1,
        help="Number of processors to use for parallel \
                        calculations.",
    )
    PARSER.add_argument(
        "--max_variants_per_compound",
        "-m",
        type=int,
        metavar="V",
        help="The maximum number of variants to create per input \
                        molecule.",
    )
    PARSER.add_argument(
        "--thoroughness",
        "-t",
        type=int,
        help="How widely to search for low-energy conformers. \
                        Larger values increase run times but can produce better \
                        results.",
    )
    PARSER.add_argument(
        "--separate_output_files",
        action="store_true",
        help="Indicates that the outputs should be split between \
                        files. If true, each output .sdf file will correspond to a \
                        single input file, but different 3D conformers will still \
                        be stored in the same file.",
    )
    PARSER.add_argument(
        "--add_pdb_output",
        action="store_true",
        help="Indicates that the outputs should also be written in \
                        the .pdb format. Creates one PDB file for each molecular \
                        variant.",
    )
    PARSER.add_argument(
        "--add_html_output",
        action="store_true",
        help="Indicates that the outputs should also be written in \
                        the .html format, for debugging. Attempts to open a \
                        browser for viewing.",
    )
    PARSER.add_argument(
        "--min_ph", metavar="MIN", type=float, help="Minimum pH to consider."
    )
    PARSER.add_argument(
        "--max_ph", metavar="MAX", type=float, help="Maximum pH to consider."
    )
    PARSER.add_argument(
        "--pka_precision",
        metavar="D",
        type=float,
        help="Size of pH substructure ranges. See Dimorphite-DL \
                        publication for details.",
    )
    PARSER.add_argument(
        "--skip_optimize_geometry",
        action="store_true",
        help="Skips the optimization step.",
    )
    PARSER.add_argument(
        "--skip_alternate_ring_conformations",
        action="store_true",
        help="Skips the non-aromatic ring-conformation \
                        generation step.",
    )
    PARSER.add_argument(
        "--skip_adding_hydrogen", action="store_true", help="Skips the ionization step."
    )
    PARSER.add_argument(
        "--skip_making_tautomers",
        action="store_true",
        help="Skips tautomer-generation step.",
    )
    PARSER.add_argument(
        "--skip_enumerate_chiral_mol",
        action="store_true",
        help="Skips the ennumeration of unspecified chiral \
                        centers.",
    )
    PARSER.add_argument(
        "--skip_enumerate_double_bonds",
        action="store_true",
        help="Skips the ennumeration of double bonds.",
    )

    PARSER.add_argument(
        "--let_tautomers_change_chirality",
        action="store_true",
        help="Allow tautomers that change \
                        the total number of chiral centers (see README.md for \
                        further explanation).",
    )

    PARSER.add_argument(
        "--use_durrant_lab_filters",
        action="store_true",
        help="Use substructure filters to \
                        remove molecular variants that, though technically \
                        possible, were judged improbable by members of the \
                        Durrant lab. See README.md for more details.",
    )

    PARSER.add_argument(
        "--2d_output_only",
        action="store_true",
        help="Skips the generate-3D-models step.",
    )
    PARSER.add_argument(
        "--cache_prerun",
        "-c",
        action="store_true",
        help="Run this before running Gypsum-DL in mpi mode.",
    )

    ARGS_DICT = vars(PARSER.parse_args())
    INPUTS = copy.deepcopy(ARGS_DICT)

    for k, v in ARGS_DICT.items():
        if v is None:
            del INPUTS[k]
    prepare_molecules(INPUTS)
    utils.log("Finished Gypsum-DL")
