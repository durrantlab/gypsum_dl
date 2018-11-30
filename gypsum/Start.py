"""
Contains the ConfGenerator object which reads, converts, and writes
small molecules.
"""

import __future__

import sys
import json
import os
from datetime import datetime
from collections import OrderedDict

import gypsum.Utils as Utils
from gypsum.parallelizer import Parallelizer

try:
    from rdkit.Chem import AllChem
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

try:
    import numpy
except:
    Utils.log("You need to install numpy and its dependencies.")
    raise ImportError("You need to install numpy and its dependencies.")

try:
    from scipy.cluster.vq import kmeans2
except:
    Utils.log("You need to install scipy and its dependencies.")
    raise ImportError("You need to install scipy and its dependencies.")

from gypsum.MolContainer import MolContainer
from gypsum.Steps.SMILES.PrepareSmiles import prepare_smiles
from gypsum.Steps.ThreeD.PrepareThreeD import prepare_three_d
from gypsum.Steps.IO.ProcessOutput import proccess_output
from gypsum.Steps.IO.LoadFiles import load_smiles_file
from gypsum.Steps.IO.LoadFiles import load_sdf_file



# see http://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
def get_user_params(args):
    """
    A function for preparing small-molecule models for docking. To work, it
    requires the python modules rdkit and molvs installed on the system.

    :param string param_file: A json file specifying the parameters.
    """
    start_time = datetime.now()

    json_warning_list = ['source', 'output_file', 'num_processors', 
                         'min_ph', 'max_ph', 'delta_ph_increment', 
                         'thoroughness', 'max_variants_per_compound',
                         'ph_std_dev']
    need_to_print_override_warning = False
    # Load the parameters from the json
    if 'json' in args:
        params = json.load(open(args['json']))
        params = set_parameters(params)
        if [i for i in json_warning_list if i in args.keys()]:
            need_to_print_override_warning = True

    else:
        params = set_parameters(args)

    if isinstance(params["source"], str):
        # smiles must be array of strs
        src = params["source"]
        if src.lower().endswith(".smi") or src.lower().endswith(".can"):
            # It's an smi file.
            smiles_data = load_smiles_file(src)
        elif params["source"].lower().endswith(".sdf"):
            # It's an sdf file. Convert it to a smiles.
            smiles_data = load_sdf_file(src)
        else:
            smiles_data = [params["source"]]
    else:
        pass  # It's already in the required format.
    params["smiles_data"] = smiles_data
    # Handle Serial overriding num_processors
    # serial fixes it to 1 processor
    if params["multithread_mode"] == "serial" or params["multithread_mode"]=="Serial":
        
        if params["num_processors"] != 1:
            print("Because --multithread_mode was set to serial, this will be run on a single processor.")
        params["num_processors"] = 1

    # Handle mpi errors if mpi4py isn't installed
    if params["multithread_mode"] == "mpi" or params["multithread_mode"] == "MPI":
        try:
            import mpi4py
        except:
            printout = "mpi4py not installed but --multithread_mode is set to mpi. \n Either install mpi4py or switch multithread_mode to multithreading or serial"
            raise ImportError(printout)

    # # # launch mpi workers
    if params["multithread_mode"] == 'mpi':
        params["Parallelizer"] = Parallelizer(params["multithread_mode"], params["num_processors"])
    else:
        # Lower level mpi (ie making a new Parallelizer within an mpi) 
        #   has problems with importing the MPI enviorment and mpi4py
        #   So we will flag it to skip the MPI mode and just go to multithread/serial
        # This is a saftey precaution
        params["Parallelizer"] = Parallelizer(params["multithread_mode"], params["num_processors"], True)

    if need_to_print_override_warning == True:
        print("WARNING: Using the --json flag overrides all other flags.")
    params["start_time"] = start_time
    return params

def conf_generator(params):
    # For Debugging
    # print("")
    # print("###########################")
    # print("num_processors  :  ", params["num_processors"])
    # print("chosen mode  :  ", params["multithread_mode"])
    # print("Parallel style:  ", params["Parallelizer"].return_mode())
    # print("Number Nodes:  ", params["Parallelizer"].return_node())
    # print("###########################")
    # print("")

    smiles_data = params["smiles_data"]

    # Make the containers
    contnrs = []
    idx_counter = 0
    for i in range(0,len(smiles_data)):
        
        smiles, name, props = smiles_data[i]
        if detect_unassigned_bonds(smiles) is None:
            print("Warning: Throwing out smile because of unassigned bonds: " + smiles)
            continue

        new_contnr = MolContainer(smiles, name, idx_counter, props)
        if new_contnr.orig_smi_canonical==None or type(new_contnr.orig_smi_canonical) !=str:
            print("Warning: Throwing out smile because of it couldn't convert to mol: " + smiles)
            continue
        
        contnrs.append(new_contnr)
        idx_counter += 1

    # Sanity check to remove None types from failed conversion
    contnrs = [x for x in contnrs if x.orig_smi_canonical!=None]
    if len(contnrs)!= idx_counter: 
        raise Exception("There is a corrupted container")

    # RUNNING MAIN SCRIPT
    prepare_smiles(contnrs, params)

    prepare_three_d(contnrs, params)

    add_mol_id_props(contnrs)
    print_current_smiles(contnrs)

    # Write any mols that fail entirely to a file.
    deal_with_failed_molecules(contnrs, params)
    start_time =params["start_time"]
    end_time = datetime.now()
    run_time = end_time - start_time
    params["start_time"] = str(start_time)
    params["end_time"] = str(end_time)
    params["run_time"] = str(run_time)

    proccess_output(contnrs, params)

    # # kill mpi workers
    params["Parallelizer"].end(params["multithread_mode"])

def detect_unassigned_bonds(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() == 0:
            return None
    return smiles

def set_parameters(params_unicode):
    """
    Set the parameters that will control this ConfGenerator object.

    :param {} params: The parameters. A dictionary of {parameter name:
                value}.
    """
    # Set the default values.
    default = OrderedDict({
        "source" : '',
        "output_file" : '',
        "separate_output_files" : False,
        "num_processors" : -1,
        "start_time" : 0,
        "end_time" : 0,
        "run_time" : 0,
        "min_ph" : 6.4,
        "max_ph" : 8.4,
        "ph_std_dev" : 1.0,
        "thoroughness" : 3,
        "max_variants_per_compound" : 5,
        "second_embed" : False,
        "2d_output_only" : False,
        "skip_optimize_geometry" : False,
        "skip_alternate_ring_conformations" : False,
        "skip_adding_hydrogen" : False,
        "skip_making_tautomers" : False,
        "skip_ennumerate_chiral_mol" : False,
        "skip_ennumerate_double_bonds" : False,
        "multithread_mode" : "multithreading"
    })

    # Modify params so that they keys are always lower case.
    # Also, rdkit doesn't play nice with unicode, so convert to ascii

    # Because Python2 & Python3 use different string objects, we separate their
    # usecases here.
    params = {}
    if sys.version_info < (3,):
        for param in params_unicode:
            val = params_unicode[param]
            if isinstance(val, unicode):
                val = str(val).encode("utf8")
            key = param.lower().encode("utf8")
            params[key] = val
    else:
        for param in params_unicode:
            val = params_unicode[param]
            key = param.lower()
            params[key] = val

    # Overlays the user parameters where they exits.
    merge_parameters(default, params)

    # Checks and prepares the final parameter list
    default = finalize_params(default)

    return default # PARAMS

def merge_parameters(default, params):
    """Combines the defaults with the user parameters."""
    # Generate a dictionary of the types
    type_dict = make_type_dict(default)

    # Move user-specified values into the parameter
    for param in params:
        # Throw an error if there's an unrecognized parameter
        if param not in default:
            Utils.log(
                "ERROR! Parameter \"" + str(param) + "\" not recognized!"
            )
            Utils.log("Here are the options:")
            Utils.log(str(default.keys()))
            raise KeyError("Unrecognized parameter")

        # Throw an error if the input parameter has a different type than
        # the default one.
        if not isinstance(params[param], type_dict[param]):
            Utils.log(
                "ERROR! The parameter \"" + param + "\" must be of " +
                "type" + str(type_dict[param]) + ", but it is of type " +
                str(type(params[param])) + "."
            )
            raise TypeError("Input parameter has a different type than the default.")

        default[param] = params[param]

def make_type_dict(dictionary):
    """Creates a dictionary of types from an existant dictionary."""
    type_dict = {}
    allowed_types = [int, float, bool, str]
    for key in dictionary:
        val = dictionary[key]
        for allowed in allowed_types:
            if isinstance(val, allowed):
                type_dict[key] = allowed
        if key not in type_dict:
            Utils.log(
                "ERROR: There appears to be an error in your parameter " +
                "JSON file. No value can have type " + str(type(val)) +
                "."
            )
            raise Exception("ERROR: There appears to be an error in your parameter")

    return type_dict

def finalize_params(dictionary):
    """Checks and updates parameters to their final values."""
    # Throw an error if there's a missing parameter.
    if dictionary["source"] == "":
        Utils.log(
            "ERROR! Missing parameter \"source\". You need to specify " +
            "the source of the input molecules (probably a SMI or SDF " +
            "file)."
        )
        raise NotImplementedError("Missing parameter")

    # Note on parameter "source", the data source. If it's a string that
    # ends in ".smi", it's treated as a smiles file. If it's a string that
    # ends in ".sdf", it's treated as an sdf file. If it's any other
    # string, it's assumed to be a smiles string itself and is assigned a
    # name of "". If it's a list, it's assumed to be a list of tuples,
    # [SMILES, Name].

    if dictionary["output_file"] == "" and dictionary["source"] != "":
        dictionary["output_file"] = dictionary["source"] + ".output.sdf"

    if dictionary["output_file"] == "":
        Utils.log(
            "ERROR! Missing parameter \"output_file\". You need to " +
            "specify where to write the output. Can be an HTML or " +
            "SDF file."
        )
        raise Exception("Missing parameter of where to write the output." +
                        "Can be an HTML or .SDF file.")

    return dictionary

def print_current_smiles(contnrs):
    """
    Prints the smiles of the current containers.
    """

    # For debugging.
    print("    Contents of MolContainers")
    for i, mol_cont in enumerate(contnrs):
        Utils.log("\t\t" + str(i) + " " + str(mol_cont.all_smiles()))
    
def add_mol_id_props(contnrs):
    """
    Once all molecules have been generated, go through each and add the
    name and a unique id (for writing to the SDF file, for example).
    """
    cont_id = 0
    for contnr in contnrs:
        for mol in contnr.mols:
            cont_id = cont_id + 1
            mol.setRDKitMolProp("UniqueID", str(cont_id))
            mol.setAllRDKitMolProps()
   
def deal_with_failed_molecules(contnrs, params):
    """
    Removes and logs failed molecules.
    """
    failed_ones = []
    for contnr in contnrs:
        if len(contnr.mols) == 0:
            astr = contnr.orig_smi + "\t" + contnr.name
            failed_ones.append(astr)

    if len(failed_ones) > 0:
        Utils.log(
            "\n3D models could not be generated for the following entries:"
        )
        Utils.log("\n".join(failed_ones))
        Utils.log("\n")

        outfile = open(params["output_file"] + ".failed.smi", 'w')
        outfile.write("\n".join(failed_ones))
        outfile.close()
