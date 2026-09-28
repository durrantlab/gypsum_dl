"""
Contains the prepare_molecules definition which reads, prepares, and writes
small molecules.
"""

from typing import Any

import json
import os
import random
import sys
from collections import OrderedDict
from datetime import datetime

import numpy
from rdkit import Chem

from gypsum_dl import utils
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.parallelizer import (
    MPI4PY_MISSING_MSG,
    MPI4PY_VERSION_MSG,
    MPI_LAUNCH_FLAG_MSG,
    Parallelizer,
    mpi4py_launch_flag_present,
    mpi4py_version_supported,
)
from gypsum_dl.steps.conf.PrepareThreeD import prepare_3d
from gypsum_dl.steps.io.LoadFiles import load_sdf_file, load_smiles_file
from gypsum_dl.steps.io.ProcessOutput import proccess_output
from gypsum_dl.steps.smiles.PrepareSmiles import prepare_smiles

# The modes Parallelizer recognizes. Anything else falls into its bare "else"
# branch, which silently substitutes multiprocessing while params keeps the
# unrecognized value, so the run either ignores the user's choice or dies at
# the first fan-out step with a message about development overrides.
VALID_JOB_MANAGERS = ("mpi", "multiprocessing", "serial")

# Ceiling on "thoroughness". The parameter multiplies into the embedding count
# in Minimize3D and GenerateAlternate3DNonaromaticRingConfs, into the
# ionization budget in AddHydrogens, and into the tautomer and double-bond
# budgets, so a mistyped value buys hours of conformer generation and the
# memory to hold it rather than an error. The default is 3 and a deliberately
# exhaustive run is still well under this, so a value above it is a typo (a
# stray digit, a shifted decimal) far more often than a request.
MAX_THOROUGHNESS = 1000


# see http://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
def prepare_molecules(args: dict[str, Any]) -> None:
    """A function for preparing small-molecule models for docking. To work, it
    requires that the python module rdkit be installed on the system.

    Args:
        args: The arguments, from the commandline.
    """

    # Keep track of the tim the program starts.
    start_time = datetime.now()

    # Whether to warn the user that the other parameters they supplied, if
    # any, will be ignored.
    need_to_print_override_warning = False

    if "json" in args:
        # "json" is one of the parameters, so we'll be ignoring the rest.
        try:
            with open(args["json"], encoding="utf-8") as json_file:
                params = json.load(json_file)
        except (OSError, ValueError):
            utils.exception("Is your input json file properly formed?")

        params = set_parameters(params)

        # The json path discards everything else the user supplied, so warn
        # about whatever is left beside "json". run.py keeps unsupplied
        # arguments out of this dictionary, so anything here was asked for:
        # previously only a hand-maintained subset of names could raise the
        # warning, and --job_manager (whose argparse default was never None)
        # could not raise it at all, so an mpi run silently became 32
        # independent multiprocessing runs writing to one output folder.
        if [key for key in args if key != "json"]:
            need_to_print_override_warning = True
    else:
        # We're actually going to use all the command-line parameters. No
        # warning necessary.
        params = set_parameters(args)

    # Seed the random number generators once, before any sampling or
    # shuffling.
    seed_random_number_generators(params)

    # If running in serial mode, make sure only one processor is used.
    if params["job_manager"] == "serial":
        if params["num_processors"] != 1:
            utils.log(
                "Because --job_manager was set to serial, this will be run on a single processor."
            )
        params["num_processors"] = 1

    # Handle mpi errors if mpi4py isn't installed
    if params["job_manager"] == "mpi":
        # The gate itself (how the launch flag is detected, which mpi4py
        # versions are acceptable, and the wording of each message) lives in
        # parallelizer, which asks the same questions when it decides whether
        # to demote mpi to multiprocessing. Only the reaction differs here:
        # the user asked for mpi explicitly, so a failed check ends the run
        # rather than quietly picking another job manager. The launch flag is
        # checked first, before mpi4py is imported, because mpi4py overrides
        # the way exceptions propagate.
        if not mpi4py_launch_flag_present():
            print(MPI_LAUNCH_FLAG_MSG)
            utils.exception(MPI_LAUNCH_FLAG_MSG)

        # Check mpi4py import
        try:
            import mpi4py
        except Exception:
            print(MPI4PY_MISSING_MSG)
            utils.exception(MPI4PY_MISSING_MSG)

        if not mpi4py_version_supported(mpi4py.__version__):
            print(MPI4PY_VERSION_MSG)
            utils.exception(MPI4PY_VERSION_MSG)

    # Throw a message if running on windows. Windows doesn't deal with with
    # multiple processors, so use only 1.
    if sys.platform == "win32":
        utils.log(
            "WARNING: Multiprocessing is not supported on Windows. Tasks will be run in Serial mode."
        )
        params["num_processors"] = 1
        params["job_manager"] = "serial"

    # Launch mpi workers if that's what's specified.
    if params["job_manager"] == "mpi":
        params["Parallelizer"] = Parallelizer(
            params["job_manager"], params["num_processors"]
        )
    else:
        # Lower-level mpi (i.e. making a new Parallelizer within an mpi) has
        # problems with importing the MPI environment and mpi4py. So we will
        # flag it to skip the MPI mode and just go to multiprocess/serial.
        # This is a saftey precaution
        params["Parallelizer"] = Parallelizer(
            params["job_manager"], params["num_processors"], True
        )

    # Let the user know that their command-line parameters will be ignored, if
    # they have specified a json file.
    if need_to_print_override_warning == True:
        utils.log("WARNING: Using the --json flag overrides all other flags.")

    # The parameters record is written during the run (by way of
    # execute_gypsum_dl), so the start time has to be in params before that
    # call rather than after it.
    params["start_time"] = str(start_time)

    # If running in mpi mode, separate_output_files must be set to true.
    if params["job_manager"] == "mpi" and params["separate_output_files"] == False:
        utils.log(
            "WARNING: Running in mpi mode, but separate_output_files is not set to True. Setting separate_output_files to True anyway."
        )
        params["separate_output_files"] = True

    # Outputing HTML files not supported in mpi mode.
    if params["job_manager"] == "mpi" and params["add_html_output"] == True:
        utils.log(
            "WARNING: Running in mpi mode, but add_html_output is set to True. HTML output is not supported in mpi mode."
        )
        params["add_html_output"] = False

    # Warn the user if he or she is not using the Durrant lab filters.
    if params["use_durrant_lab_filters"] is False:
        utils.log(
            "WARNING: Running Gypsum-DL without the Durrant-lab filters. In looking over many Gypsum-DL-generated "
            + "variants, we have identified a number of substructures that, though technically possible, strike us "
            + "as improbable or otherwise poorly suited for virtual screening. We strongly recommend removing these "
            + "by running Gypsum-DL with the --use_durrant_lab_filters option.",
            trailing_whitespace="\n",
        )

    # Load SMILES data. finalize_params has already established that "source"
    # is a string naming an existing file, so only the extension is left to
    # dispatch on.
    utils.log("Loading molecules from " + os.path.basename(params["source"]) + "...")

    src = params["source"]
    if src.lower().endswith(".smi") or src.lower().endswith(".can"):
        # It's an smi file.
        smiles_data = load_smiles_file(src)
    elif src.lower().endswith(".sdf"):
        # It's an sdf file. Convert it to a smiles.
        smiles_data = load_sdf_file(src)
    else:
        # Raising here keeps the diagnosis on the actual problem. Falling
        # through instead left smiles_data holding the filename itself, which
        # only failed later, while unpacking a (smiles, name, props) tuple.
        utils.exception(
            'Your "source" parameter must name a file that ends in a .can, '
            f".smi, or .sdf extension, but it names {os.path.basename(src)}."
        )

    # Make the output directory if necessary.
    try:
        os.makedirs(params["output_folder"], exist_ok=True)
    except OSError:
        utils.exception("Output folder directory couldn't be found or created.")

    # For Debugging
    # print("")
    # print("###########################")
    # print("num_procs  :  ", params["num_processors"])
    # print("chosen mode  :  ", params["job_manager"])
    # print("Parallel style:  ", params["Parallelizer"].return_mode())
    # print("Number Nodes:  ", params["Parallelizer"].return_node())
    # print("###########################")
    # print("")

    # Make the molecule containers. The index is the position of the record in
    # the input file, not a count of the records accepted: every external
    # number (output filenames, UniqueID, the failed-molecule file) derives
    # from it, so it has to line up with the input the user can look at. A
    # rejected record therefore leaves a gap, which nothing downstream minds
    # (utils.contnrs_by_idx needs only uniqueness, and mpi mode rewrites the
    # live index to 0 anyway).
    contnrs = []
    for i in range(0, len(smiles_data)):
        try:
            smiles, name, props = smiles_data[i]
        except Exception:
            msg = 'Unexpected error. Does your "source" parameter specify a '
            msg = msg + "filename that ends in a .can, .smi, or .sdf extension?"
            utils.exception(msg)

        if detect_unassigned_bonds(smiles) is None:
            utils.log(
                "WARNING: Throwing out SMILES because of unassigned bonds: " + smiles
            )
            continue

        new_contnr = MolContainer(smiles, name, i, props)
        if (
            new_contnr.orig_smi_canonical == None
            or type(new_contnr.orig_smi_canonical) != str
        ):
            utils.log(
                "WARNING: Throwing out SMILES because of it couldn't convert to mol: "
                + smiles
            )
            continue

        contnrs.append(new_contnr)

    # Remove None types from failed conversion
    contnrs = [x for x in contnrs if x.orig_smi_canonical != None]

    # In multiprocessing mode, Gypsum-DL parallelizes each small-molecule
    # preparation step separately. But this scheme is inefficient in MPI mode
    # because it increases the amount of communication required between nodes.
    # So for MPI mode, we will run all the preparation steps for a given
    # molecule container on a single thread.
    if params["Parallelizer"].return_mode() != "mpi":
        # Non-MPI (e.g., multiprocessing)
        execute_gypsum_dl(contnrs, params)
    else:
        # MPI mode. Group the molecule containers so they can be passed to the
        # parallelizer.
        job_input = []
        temp_param = {}
        for key in list(params.keys()):
            if key == "Parallelizer":
                temp_param["Parallelizer"] = None
            else:
                temp_param[key] = params[key]

        for contnr in contnrs:
            # Each container is run in isolation, so it becomes container zero
            # of its own job. Go through update_idx rather than assigning the
            # attribute: mol_orig_frm_inp_smi must be restamped too, because
            # the desalter hands that very object to the container for
            # single-fragment inputs, and later steps regroup molecules by
            # mol.contnr_idx.
            contnr.update_idx(0)
            job_input.append(tuple([[contnr], temp_param]))
        job_input = tuple(job_input)

        params["Parallelizer"].run(job_input, execute_gypsum_dl)

    # Calculate the total run time. Neither of these is a parameter: the
    # parameters record is written by save_to_sdf while the run is still going,
    # so it can carry start_time but never these two. They are reported here
    # only.
    end_time = datetime.now()
    run_time = end_time - start_time

    utils.log("\nStart time at: " + str(start_time))
    utils.log("End time at:   " + str(end_time))
    utils.log("Total time at: " + str(run_time))

    # Kill mpi workers if necessary.
    params["Parallelizer"].end(params["job_manager"])


def seed_random_number_generators(params: dict[str, Any]) -> None:
    """Seed both of the global random number generators Gypsum-DL draws from.

    Variant selection samples with the random module, and the ring-conformation
    clustering samples with numpy (by way of scipy's kmeans2), so seeding one
    without the other still leaves a seeded run varying between invocations.
    Seeding reaches only the calling process, which makes serial runs
    reproducible; multiprocessing and mpi runs also depend on how tasks land on
    workers and on each worker's own generator state, so they remain
    nondeterministic.

    Args:
        params: The parameters, which may carry a non-negative random_seed. A
            negative seed (the default) leaves both generators alone.
    """

    seed = params.get("random_seed", -1)
    if seed < 0:
        return

    random.seed(seed)

    # numpy rejects seeds that do not fit in 32 bits, while the random module
    # accepts an integer of any size. Fold larger values in rather than
    # failing the run over the choice of seed.
    numpy.random.seed(seed % 2**32)

    utils.log(
        "Using random_seed = "
        + str(seed)
        + ". Note that this makes serial runs reproducible; multiprocessing "
        + "and mpi runs remain nondeterministic."
    )


def execute_gypsum_dl(contnrs: list, params: dict[str, Any]) -> None:
    """A function for doing all of the manipulations to each molecule.

    Args:
        contnrs: A list of all molecules.
        params: A dictionary containing all of the parameters.
    """
    # Start creating the models.

    # Prepare the smiles. Desalt, consider alternate ionization, tautometeric,
    # stereoisomeric forms, etc.
    prepare_smiles(contnrs, params)

    # Convert the processed SMILES strings to 3D.
    prepare_3d(contnrs, params)

    # Add in name and unique id to each molecule.
    add_mol_id_props(contnrs)

    # Output the current SMILES. This is a diagnostic dump of every variant of
    # every container, and it forces a hydrogen-stripped canonicalization per
    # variant purely to print it, so it follows the same debug flag as the
    # identical calls in prepare_smiles.
    if params.get("debug", False):
        utils.print_current_smiles(contnrs)

    # Write any mols that fail entirely to a file.
    deal_with_failed_molecules(contnrs, params)  ####

    # Process the output.
    proccess_output(contnrs, params)


def detect_unassigned_bonds(smiles: str) -> str | None:
    """Detects whether a give smiles string has unassigned bonds.

    Args:
        smiles: The smiles string.

    Returns:
        None if it has bad bonds, or the input smiles string otherwise.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        # Apparently the bonds are particularly bad, because couldn't even
        # create the molecule.
        return None
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 0:
            return None
    return smiles


def set_parameters(params_unicode: dict[str, Any]) -> dict[str, Any]:
    """Set the parameters that will control this ConfGenerator object.

    Args:
        params_unicode: The parameters, with keys and values possibly in
            unicode.

    Returns:
        The parameters, properly processed, with defaults used when no
            value specified.
    """

    # Set the default values.
    default = OrderedDict(
        {
            "source": "",
            "output_folder": "./",
            "separate_output_files": False,
            "add_pdb_output": False,
            "add_html_output": False,
            "num_processors": -1,
            "start_time": 0,
            "min_ph": 6.4,
            "max_ph": 8.4,
            "pka_precision": 1.0,
            "thoroughness": 3,
            "max_variants_per_compound": 5,
            "second_embed": False,
            "2d_output_only": False,
            "skip_optimize_geometry": False,
            "skip_alternate_ring_conformations": False,
            "skip_adding_hydrogen": False,
            "skip_making_tautomers": False,
            "skip_enumerate_chiral_mol": False,
            "skip_enumerate_double_bonds": False,
            "let_tautomers_change_chirality": False,
            "use_durrant_lab_filters": False,
            "job_manager": "multiprocessing",
            "test": False,
            # Gates the diagnostic container dumps in prepare_smiles and
            # execute_gypsum_dl. It has to be listed here, because
            # merge_parameters rejects any key missing from the defaults.
            "debug": False,
            # Seed for the global random and numpy generators. A value >= 0
            # makes serial runs reproducible; multiprocessing and mpi runs
            # remain nondeterministic. A negative value leaves both
            # generators unseeded (previous behavior).
            "random_seed": -1,
        }
    )

    # Modify params so that they keys are always lower case.
    # Also, rdkit doesn't play nice with unicode, so convert to ascii

    # Because Python2 & Python3 use different string objects, we separate their
    # usecases here.
    params = {}
    for param in params_unicode:
        val = params_unicode[param]
        key = param.lower()
        params[key] = val

    # Overwrites values with the user parameters where they exit.
    merge_parameters(default, params)

    # Checks and prepares the final parameter list.
    return finalize_params(default)


def merge_parameters(default: dict[str, Any], params: dict[str, Any]) -> None:
    """Add default values if missing from parameters.

    Args:
        default: The parameters.
        params: The default values
    """

    # Generate a dictionary with the same keys, but the types for the values.
    type_dict = make_type_dict(default)

    # Move user-specified values into the parameter.
    for param in params:
        # Throw an error if there's an unrecognized parameter.
        if param not in default:
            utils.log(f'Parameter "{str(param)}" not recognized!')
            utils.log("Here are the options:")
            utils.log(" ".join(sorted(list(default.keys()))))
            utils.exception(f"Unrecognized parameter: {str(param)}")

        # Throw an error if the input parameter has a different type than
        # the default one.
        if not isinstance(params[param], type_dict[param]):
            # Cast int to float if necessary
            if type(params[param]) is int and type_dict[param] is float:
                params[param] = float(params[param])
            else:
                # Seems to be a type mismatch.
                utils.exception(
                    'The parameter "'
                    + param
                    + '" must be of '
                    + "type "
                    + str(type_dict[param])
                    + ", but it is of type "
                    + str(type(params[param]))
                    + "."
                )

        # Update the parameter value with the user-defined one.
        default[param] = params[param]


def make_type_dict(dictionary: dict[str, Any]) -> dict[str, Any]:
    """Creates a types dictionary from an existant dictionary. Keys are
    preserved, but values are the types.

    Args:
        dictionary: A dictionary, with keys are values.

    Returns:
        A dictionary with the same keys, but the values are the types.
    """
    type_dict = {}
    allowed_types = [int, float, bool, str]
    # Go through the dictionary keys.
    for key in dictionary:
        # Get the the type of the value.
        val = dictionary[key]
        for allowed in allowed_types:
            if isinstance(val, allowed):
                # Add it to the type_dict.
                type_dict[key] = allowed

        # The value ha san unacceptable type. Throw an error.
        if key not in type_dict:
            utils.exception(
                "ERROR: There appears to be an error in your parameter "
                + "JSON file. No value can have type "
                + str(type(val))
                + "."
            )

    return type_dict


def finalize_params(params: dict[str, Any]) -> dict[str, Any]:
    """Checks and updates parameters to their final values.

    Args:
        params: The parameters.

    Returns:
        The parameters, corrected/updated where needed.
    """

    # Throw an error if there's a missing parameter.
    if params["source"] == "":
        utils.exception(
            'Missing parameter "source". You need to specify '
            + "the source of the input molecules (probably a SMI or SDF "
            + "file)."
        )

    # Note on parameter "source", the data source. It must name a file: one
    # ending in ".smi" or ".can" is treated as a SMILES file, and one ending in
    # ".sdf" as an SDF file.

    # Check some required variables. Note that os.path.abspath does not touch
    # the filesystem and so never reports a missing file; the existence check
    # has to be made separately.
    params["source"] = os.path.abspath(params["source"])
    if not os.path.isfile(params["source"]):
        utils.exception(f"Source file not found: {params['source']}")
    source_dir = os.path.dirname(params["source"]) + os.sep

    # An empty output_folder is always filled in here, either from the source
    # directory or from the "./" default applied when the parameters are
    # merged, so the .pdb and separate-file output modes cannot reach this
    # point without a folder to write to.
    if params["output_folder"] == "" and params["source"] != "":
        params["output_folder"] = f"{source_dir}output{str(os.sep)}"

    # if not os.path.exists(params["output_folder"]) or not os.path.isdir(params["output_folder"]):
    #     utils.exception(
    #         "The specified \"output_folder\", " + params["output_folder"] +
    #         ", either does not exist or is a file rather than a folder. " +
    #         "Please provide the path to an existing folder instead."
    #     )

    # Range-check the two parameters that control how many candidates each
    # step generates. Only their types are checked when they are merged, and
    # out-of-range values fail far from here: thoroughness reaches
    # math.log(thoroughness * max_variants_per_compound, 2) inside a worker,
    # where a zero raises and the failure is swallowed as a dropped molecule,
    # and a negative value of either reaches list slices that quietly truncate
    # the candidate pools instead of raising. Guard each check with a
    # membership test, because finalize_params is also called with partial
    # parameter dictionaries.
    if "thoroughness" in params and params["thoroughness"] < 1:
        utils.exception('The parameter "thoroughness" must be at least 1.')

    if "thoroughness" in params and params["thoroughness"] > MAX_THOROUGHNESS:
        utils.exception(
            'The parameter "thoroughness" must be '
            + str(MAX_THOROUGHNESS)
            + " or less, but it is "
            + str(params["thoroughness"])
            + "."
        )

    # Zero is allowed for max_variants_per_compound: the SMILES enumeration
    # steps treat it as a sentinel that skips them entirely.
    if (
        "max_variants_per_compound" in params
        and params["max_variants_per_compound"] < 0
    ):
        utils.exception(
            'The parameter "max_variants_per_compound" must be 0 or greater.'
        )

    # Make sure job_manager is always lower case.
    params["job_manager"] = params["job_manager"].lower()

    # Only the command-line path restricts this parameter (through argparse
    # choices); the json and API paths reach here with any string at all.
    if params["job_manager"] not in VALID_JOB_MANAGERS:
        utils.exception(
            'The parameter "job_manager" must be one of '
            + ", ".join(VALID_JOB_MANAGERS)
            + ', but it is "'
            + params["job_manager"]
            + '".'
        )

    # An inverted pH window is passed straight to Dimorphite-DL, which has no
    # reason to expect one. Guard with a membership test, because
    # finalize_params is also called with partial parameter dictionaries.
    if (
        "min_ph" in params
        and "max_ph" in params
        and params["min_ph"] > params["max_ph"]
    ):
        utils.exception(
            'The parameter "min_ph" ('
            + str(params["min_ph"])
            + ') cannot be greater than "max_ph" ('
            + str(params["max_ph"])
            + ")."
        )

    return params


def add_mol_id_props(contnrs: list[MolContainer]) -> None:
    """Once all molecules have been generated, go through each and add the
       name and a unique id (for writing to the SDF file, for example).

    Args:
        contnrs: A list of containers (MolContainer.MolContainer).
    """

    # A plain running counter is only unique within one call, and every mpi
    # rank calls this independently (as does every task in
    # separate_output_files mode), so concatenating the per-input files gave
    # many records the same id. Qualifying the per-container variant number
    # with the container's original input index makes the id unique across the
    # whole run, and matches how the separate output files are named.
    # The id goes into mol_props rather than straight onto the rdkit mol so
    # that add_container_properties, which fills gaps in mol_props from the
    # input file, cannot substitute an input UniqueID tag for the id the PDB
    # and SDF filenames are built to match.
    for contnr in contnrs:
        for variant_id, mol in enumerate(contnr.mols, start=1):
            mol.mol_props["UniqueID"] = f"{contnr.contnr_idx_orig + 1}_{variant_id}"
            mol.set_all_rdkit_mol_props()


def deal_with_failed_molecules(
    contnrs: list[MolContainer], params: dict[str, Any]
) -> None:
    """Removes and logs failed molecules.

    Args:
        contnrs: A list of containers (MolContainer.MolContainer).
        params: The parameters, used to determine the filename that will
            contain the failed molecules.
    """

    # To keep track of failed molecules
    failed = [contnr for contnr in contnrs if len(contnr.mols) == 0]
    if not failed:
        return

    lines: dict[MolContainer, str] = {
        contnr: f"{contnr.orig_smi}\t{contnr.name}" for contnr in failed
    }

    # Let the user know if there's more than one failed molecule.
    utils.log("\n3D models could not be generated for the following entries:")
    utils.log("\n".join(lines.values()))
    utils.log("\n")

    # Write the failures to an smi file. In separate-file mode (which mpi mode
    # forces) every task would otherwise truncate and rewrite one shared path;
    # under mpi the ranks do so concurrently, leaving only the last writer's
    # failures behind. Qualifying the filename with the container's original
    # index matches the convention used for the other separate output files,
    # but it has to be done per container: an mpi task holds a single
    # container, while a non-mpi run passes all of them in one call, so keying
    # off contnrs[0] alone filed every failure under input 1.
    groups: dict[str, list[MolContainer]] = {}
    if params.get("separate_output_files", False):
        for contnr in failed:
            failed_name = f"gypsum_dl_failed__input{contnr.contnr_idx_orig + 1}.smi"
            groups.setdefault(failed_name, []).append(contnr)
    else:
        groups["gypsum_dl_failed.smi"] = failed

    for failed_name, group in groups.items():
        with open(
            os.path.join(params["output_folder"], failed_name),
            "w",
            encoding="utf-8",
        ) as outfile:
            outfile.write("\n".join(lines[contnr] for contnr in group))
