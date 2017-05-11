"""
Contains the ConfGenerator object which reads, converts, and writes
small molecules.
"""
import sys
import json
import os
from collections import OrderedDict

from gypsum import Utils

try:
    from rdkit.Chem import AllChem
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)

try:
    import numpy
except:
    Utils.log("You need to install numpy and its dependencies.")
    sys.exit(0)

try:
    from scipy.cluster.vq import kmeans2
except:
    Utils.log("You need to install scipy and its dependencies.")
    sys.exit(0)

from gypsum.MolContainer import MolContainer
from gypsum.Steps.SMILES import prepare_smiles
from gypsum.Steps.ThreeD import prepare_three_d
from gypsum.Steps.IO import load_smiles_file, load_sdf_file, proccess_output

# see http://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
class ConfGenerator(object):
    """
    A class for preparing small-molecule models for docking. To work, it
    requires the python modules rdkit and molvs, as well as openbabel
    installed as an executable on the system.
    """
    def __init__(self, args):
        """
        The class constructor.

        :param string param_file: A json file specifying the parameters.
        """
        warning_list = ['source', 'output_file', 'openbabel_executable',
                        'num_processors', 'min_ph', 'max_ph',
                        'delta_ph_increment', 'thoroughness',
                        'max_variants_per_compound']

        # Load the parameters from the json
        if args.has_key('json'):
            params = json.load(open(args['json']))
            self.set_parameters(params)
            if [i for i in warning_list if i in args.keys()]:
                print "WARNING: Using the --json flag overrides all other flags."
        else:
            self.set_parameters(args)

        if isinstance(self.params["source"], basestring):
            # smiles must be array of strs
            src = self.params["source"]
            if src.lower().endswith(".smi") or src.lower().endswith(".can"):
                # It's an smi file.
                smiles_data = load_smiles_file(src)
            elif self.params["source"].lower().endswith(".sdf"):
                # It's an sdf file. Convert it to a smiles.
                smiles_data = load_sdf_file(src)
            else:
                smiles_data = [self.params["source"]]
        else:
            pass  # It's already in the required format.

        # Make the containers
        self.contnrs = []
        for idx, data in enumerate(smiles_data):
            smiles, name = data
            new_contnr = MolContainer(smiles, name, idx)
            self.contnrs.append(new_contnr)

    def run(self):
        """
        Starts the conversion process, ultimately writing the converted files to disk.
        """
        prepare_smiles(self)

        prepare_three_d(self)

        self.add_mol_id_props()
        self.print_current_smiles()

        # Write any mols that fail entirely to a file.
        self.deal_with_failed_molecules()

        proccess_output(self)


    def set_parameters(self, params_unicode):
        """
        Set the parameters that will control this ConfGenerator object.

        :param {} params: The parameters. A dictionary of {parameter name:
                  value}.
        """
        # Set the default values.
        default = OrderedDict({
            "source" : "",
            "output_file" : "",
            "separate_output_files" : False,
            "openbabel_executable" : "/usr/local/bin/obabel",
            "num_processors" : 1,
            "min_ph" : 5.0,
            "max_ph" : 9.0,
            "delta_ph_increment" : 0.5,
            "thoroughness" : 3,
            "max_variants_per_compound" : 5,
            "skip_optimize_geometry" : False,
            "skip_alternate_ring_conformations" : False,
            "skip_adding_hydrogen" : False,
            "skip_making_tautomers" : False,
            "skip_ennumerate_chiral_mol" : False,
            "skip_ennumerate_double_bonds" : False,
            "2d_output_only" : False
        })

        # Modify params so that they keys are always lower case.
        # Also, rdkit doesn't play nice with unicode, so convert to ascii
        params = {}
        for param in params_unicode:
            val = params_unicode[param]
            if isinstance(val, basestring):
                val = val.encode("utf8")
            key = param.lower().encode("utf8")
            params[key] = val

        # Overlays the user parameters where they exits.
        self.merge_parameters(default, params)

        # Checks and prepares the final parameter list
        default = self.finalize_params(default)

        self.params = default

    def merge_parameters(self, default, params):
        """Combines the defaults with the user parameters."""
        # Generate a dictionary of the types
        type_dict = self.make_type_dict(default)

        # Move user-specified values into the parameter
        for param in params:
            # Throw an error if there's an unrecognized parameter
            if param not in default:
                Utils.log(
                    "ERROR! Parameter \"" + param + "\" not recognized!"
                )
                Utils.log("Here are the options:")
                Utils.log(str(default.keys()))
                sys.exit(0)

            # Throw an error if the input parameter has a different type that
            # the default one.
            if not isinstance(params[param], type_dict[param]):
                Utils.log(
                    "ERROR! The parameter \"" + param + "\" must be of " +
                    "type" + str(type_dict[param]) + ", but it is of type " +
                    str(type(params[param])) + "."
                )
                sys.exit(0)

            default[param] = params[param]

    @staticmethod
    def make_type_dict(dictionary):
        """Creates a dictionary of types from an existant dictionary."""
        type_dict = {}
        allowed_types = [int, float, str, bool]
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
                sys.exit(0)

        return type_dict

    @staticmethod
    def finalize_params(dictionary):
        """Checks and updates parameters to their final values."""
        # Throw an error if there's a missing parameter.
        if dictionary["source"] == "":
            Utils.log(
                "ERROR! Missing parameter \"source\". You need to specify " +
                "the source of the input molecules (probably a SMI or SDF " +
                "file)."
            )
            sys.exit(0)

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
            sys.exit(0)

        if not os.path.exists(dictionary["openbabel_executable"]):
            Utils.log(
                "ERROR! There is no executable at " +
                dictionary["openbabel_executable"] + ". Please specify the " +
                "correct path in your parameters file and/or install Open " +
                "Babel if necessary."
            )
            sys.exit(0)

        return dictionary

    def print_current_smiles(self):
        """
        Prints the smiles of the current containers.
        """

        # For debugging.
        print "    Contents of MolContainers"
        for i, mol_cont in enumerate(self.contnrs):
            Utils.log("\t\t" + str(i) + " " + str(mol_cont.all_smiles()))


    def add_mol_id_props(self):
        """
        Once all molecules have been generated, go through each and add the
        name and a unique id (for writing to the SDF file, for example).
        """
        cont_id = 0
        for contnr in self.contnrs:
            for mol in contnr.mols:
                cont_id = cont_id + 1
                mol.setRDKitMolProp("UniqueID", str(cont_id))
                mol.setAllRDKitMolProps()

    def add_indexed_mols_to_mols(self, items):
        """
        Adds a molecule to the specified MolContainer.

        :param list items: A list of tuples, [(index, mol), (index, mol), ...]
        """
        for index, mol in items:
            self.contnrs[index].add_mol(mol)


    def deal_with_failed_molecules(self):
        """
        Removes and logs failed molecules.
        """
        failed_ones = []
        for contnr in self.contnrs:
            if len(contnr.mols) == 0:
                astr = contnr.orig_smi + "\t" + contnr.name
                failed_ones.append(astr)

        if len(failed_ones) > 0:
            Utils.log(
                "\n3D models could not be generated for the following entries:"
            )
            Utils.log("\n".join(failed_ones))
            Utils.log("\n")

            outfile = open(self.params["output_file"] + ".failed.smi", 'w')
            outfile.write("\n".join(failed_ones))
            outfile.close()
