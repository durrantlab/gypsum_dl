import Multiprocess as mp
import Utils
import ChemUtils
from MolContainer import MolContainer
import Steps
from Steps.IO import load_smiles_file
from Steps.IO import load_sdf_file
import sys
import json
from collections import OrderedDict
import MyMol
import os

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


# see http://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
class ConfGenerator:
    """
    A class for preparing small-molecule models for docking. To work, it
    requires the python modules rdkit and molvs, as well as openbabel
    installed as an executable on the system.
    """

    def __init__(self, args):
        """
        The class constructor. Starts the conversion process, ultimately
        writing the converted files to disk.

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
                smiles_data = load_smiles_file(self.params["source"])
            elif self.params["source"].lower().endswith(".sdf"):
                # It's an sdf file. Convert it to a smiles.
                smiles_data = load_sdf_file(self.params["source"])
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

        # ESSENTIAL
        Steps.SMILES.desalt_orig_smi(self)

        if self.params["skip_adding_hydrogen"] == False:
            Steps.SMILES.add_hydrogens(self)
        else:
            # PUTTING STUFF HERE THAT PATRICK WILL MOVE INTO DEF OR EXTERNAL
            # MODULE PER HIS WISDOM

            # Problem: Each molecule container holds one smiles string
            # (corresponding to the input structure). obabel produces multiple
            # smiles strings at different pH values in the previous step.
            # There is no way to store muliple smiles in a molecule container.
            # But those containers are designed to store multiple RDKit
            # molecule objects. To the previous step stores the differently
            # protonated models as those objects, in the container's mol list.

            # But, if the user skips the previous step, then the one smiles
            # needs to be converted to a RDKit mol object for subsequent steps
            # to work. Let's do that here.
            
            for i, mol_cont in enumerate(self.contnrs):
                if len(mol_cont.mols) == 0:
                    smi = mol_cont.orig_smi_canonical
                    mol_cont.add_smiles(smi)

        self.print_current_smiles()

        # Do tautomers first, because obliterates chiral info I think
        if self.params["skip_making_tautomers"] == False:
            Steps.SMILES.make_tauts(self)
        self.print_current_smiles()

        if self.params["skip_ennumerate_chiral_mol"] == False:
            Steps.SMILES.enumerate_chiral_molecules(self)
        self.print_current_smiles()

        # Suprized you have a hard time generating enantiomers here:
        # CCC(C)NC(=O)CC(C)C

        if self.params["skip_ennumerate_double_bonds"] == False:
            Steps.SMILES.enumerate_double_bonds(self)
        self.print_current_smiles()

        if self.params["2d_output_only"] == False:
            Steps.ThreeD.convert_2d_to_3d(self)
        self.print_current_smiles()

        if self.params["skip_alternate_ring_conformations"] == False:
            Steps.ThreeD.generate_alternate_3d_nonaromatic_ring_confs(self)
        self.print_current_smiles()

        if self.params["skip_optimize_geometry"] == False:
            Steps.ThreeD.minimize_3d(self)
        self.print_current_smiles()
        
        self.add_mol_id_props()
        self.print_current_smiles()

        # Write any mols that fail entirely to a file.
        self.deal_with_failed_molecules()

        if self.params["output_file"].lower().endswith(".html"):
            Steps.IO.web_2d_output(self)
        else:
            Steps.IO.save_to_sdf(self)


    def set_parameters(self, params):
        """
        Set the parameters that will control this ConfGenerator object.

        :param {} params: The parameters. A dictionary of {parameter name:
                  value}.
        """

        # Set the default values.
        default = OrderedDict({})
        default["source"] = ""
        default["output_file"] = ""
        default["separate_output_files"] = False
        default["openbabel_executable"] = "/usr/local/bin/obabel"
        default["num_processors"] = 1

        default["min_ph"] = 5.0
        default["max_ph"] = 9.0
        default["delta_ph_increment"] = 0.5

        default["thoroughness"] = 3
        default["max_variants_per_compound"] = 5

        default["skip_optimize_geometry"] = False
        default["skip_alternate_ring_conformations"] = False
        default["skip_adding_hydrogen"] = False
        default["skip_making_tautomers"] = False
        default["skip_ennumerate_chiral_mol"] = False
        default["skip_ennumerate_double_bonds"] = False

        default["2d_output_only"] = False

        # Modify params so that they keys are always lower case.
        # Also, rdkit doesn't play nice with unicode, so convert to ascii
        params2 = {}
        for param in params:
            val = params[param]
            
            if isinstance(val, basestring):
                val = val.encode("utf8")
            
            params2[param.lower().encode("utf8")] = val

        # Generate a dictionary of the types
        type_dict = {}
        for key in default:
            val = default[key]
            if type(val) is int:
                type_dict[key] = int
            elif type(val) is float:
                type_dict[key] = float
            elif type(val) is str:
                type_dict[key] = str
            elif type(val) is bool:
                type_dict[key] = bool
            else:
                Utils.log(
                    "ERROR: There appears to be an error in your parameter " +
                    "JSON file. No value can have type " + str(type(val)) + 
                    "."
                )
                sys.exit(0)
        
        # Move user-specified values into the parameter
        for param in params2:
            param = param.lower()

            # Throw an error if there's an unrecognized parameter
            if not (param in default):
                Utils.log(
                    "ERROR! Parameter \"" + param + "\" not recognized!"
                )
                Utils.log("Here are the options:")
                Utils.log(str(default.keys()))
                sys.exit(0)
            
            # Throw an error if the input parameter has a different type that
            # the default one.
            if type(params2[param]) is not type_dict[param]:
                Utils.log(
                    "ERROR! The parameter \"" + param + "\" must be of " +
                    "type" + str(type_dict[param]) + ", but it is of type " + 
                    str(type(params2[param])) + "."
                )
                sys.exit(0)
            
            default[param] = params2[param]

        # Throw an error if there's a missing parameter.
        if default["source"] == "":
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

        if default["output_file"] == "" and default["source"] != "":
            default["output_file"] = default["source"] + ".output.sdf"

        if default["output_file"] == "":
            Utils.log(
                "ERROR! Missing parameter \"output_file\". You need to " +
                "specify where to write the output. Can be an HTML or " +
                "SDF file."
            )
            sys.exit(0)

        if not os.path.exists(default["openbabel_executable"]):
            Utils.log(
                "ERROR! There is no executable at " +
                default["openbabel_executable"] + ". Please specify the " +
                "correct path in your parameters file and/or install Open " +
                "Babel if necessary."
            )
            sys.exit(0)

        self.params = default


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

        id = 0
        for contnr_idx, contnr in enumerate(self.contnrs):
            for mol_index, mol in enumerate(contnr.mols):
                id = id + 1
                mol.setRDKitMolProp("UniqueID", str(id))
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
        What does this function do?
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

            f = open(self.params["output_file"] + ".failed.smi", 'w')
            f.write("\n".join(failed_ones))
            f.close()

