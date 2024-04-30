"""
A module for loading in files.
"""

from . import utils as Utils

try:
    from rdkit import Chem
except Exception:
    Utils.exception("You need to install rdkit and its dependencies.")


def load_smiles_file(filename):
    """Loads a smiles file.

    :param filename: The filename.
    :type filename: str
    :return: A list of tuples, (SMILES, Name).
    :rtype: list
    """

    # A smiles file contains one molecule on each line. Each line is a string,
    # separated by white space, followed by the molecule name.
    data = []
    duplicate_names = {}
    line_counter = 0
    name_list = []
    for line in open(filename):
        # You've got the line.
        line = line.strip()
        if line != "":
            # From that line, get the smiles string and name.
            chunks = line.split()
            smiles = chunks[0]
            name = " ".join(chunks[1:])

            # Handle unnamed ligands.
            if not name:
                name = f"untitled_line_{line_counter + 1}"
                Utils.log(
                    (
                        "\tUntitled ligand on line {}. Naming that ligand "
                        + "{}. All associated files will be referred to with "
                        + "this name."
                    ).format(line_counter + 1, name)
                )

            # Handle duplicate ligands in same list.
            if name in name_list:
                # If multiple names...
                if name in list(duplicate_names.keys()):
                    duplicate_names[name] = duplicate_names[name] + 1
                else:
                    duplicate_names[name] = 2
                new_name = f"{name}_copy_{duplicate_names[name]}"
                Utils.log(f"\nMultiple entries with the ligand name: {name}")
                Utils.log(
                    f"\tThe version of the ligand on line {line_counter} will be retitled {new_name}"
                )
                Utils.log("\tAll associated files will be referred to with this name")
                name = new_name
            # Save the data for this line and advance.
            name_list.append(name)
            line_counter += 1
            data.append((smiles, name, {}))

    # Return the data.
    return data


def load_sdf_file(filename):
    """Loads an sdf file.

    :param filename: The filename.
    :type filename: str
    :return: A list of tuples, (SMILES, Name).
    :rtype: list
    """

    suppl = Chem.SDMolSupplier(filename)
    data = []
    duplicate_names = {}
    missing_name_counter = 0
    mol_obj_counter = 0
    name_list = []
    for mol in suppl:
        # Convert mols to smiles. That's what the rest of the program is
        # designed to deal with.
        if mol:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        else:
            Utils.log(
                "\tWarning: Could not convert some SDF-formatted files to SMILES. Consider using an SMI (SMILES) file instead."
            )
            continue

        try:
            name = mol.GetProp("_Name")
        except Exception:
            name = ""

        # Handle unnamed ligands
        if not name:
            Utils.log(
                f"\tUntitled ligand for the {mol_obj_counter} molecule in the input SDF"
            )
            name = f"untitled_{missing_name_counter}_molnum_{mol_obj_counter}"
            Utils.log(f"\tNaming that ligand {name}")
            Utils.log("\tAll associated files will be referred to with this name")
            missing_name_counter += 1

            # Handle duplicate ligands in same list.
            if name in name_list:
                Utils.log(f"\nMultiple entries with the ligand name: {name}")
                # If multiple names.
                if name in list(duplicate_names.keys()):
                    duplicate_names[name] = duplicate_names[name] + 1

                else:
                    duplicate_names[name] = 2
                new_name = f"{name}_copy_{duplicate_names[name]}"
                name = new_name
                Utils.log(
                    f"\tThe version of the ligand for the {mol_obj_counter} molecule in the SDF file will be retitled {name}"
                )
                Utils.log("\tAll associated files will be referred to with this name")
            mol_obj_counter += 1
            name_list.append(name)

        # SDF files may also contain properties. Get those as well.
        try:
            properties = mol.GetPropsAsDict()
        except Exception:
            properties = {}

        if smiles != "":
            data.append((smiles, name, properties))

    return data
