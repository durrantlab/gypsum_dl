"""
A module for loading in files.
"""

from gypsum_dl import utils

try:
    from rdkit import Chem
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")


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
    name_set: set[str] = set()
    with open(filename, encoding="utf-8") as f:
        # enumerate over the raw file so reported line numbers match the file,
        # regardless of any blank lines that get skipped.
        for line_num, line in enumerate(f, start=1):
            # You've got the line.
            line = line.strip()
            if line == "":
                continue

            # From that line, get the smiles string and name.
            chunks = line.split()
            smiles = chunks[0]
            name = " ".join(chunks[1:])

            # Handle unnamed ligands.
            if not name:
                name = f"untitled_line_{line_num}"
                utils.log(
                    (
                        "\tUntitled ligand on line {}. Naming that ligand "
                        + "{}. All associated files will be referred to with "
                        + "this name."
                    ).format(line_num, name)
                )

            # Handle duplicate ligands in same list.
            if name in name_set:
                # Keep incrementing the counter until the generated name is one
                # nothing else has claimed. An input that already contains,
                # say, "lig_copy_2" would otherwise collide with the name
                # generated for a duplicate "lig".
                duplicate_names[name] = duplicate_names.get(name, 1) + 1
                new_name = f"{name}_copy_{duplicate_names[name]}"
                while new_name in name_set:
                    duplicate_names[name] += 1
                    new_name = f"{name}_copy_{duplicate_names[name]}"
                utils.log(f"\nMultiple entries with the ligand name: {name}")
                utils.log(
                    f"\tThe version of the ligand on line {line_num} will be retitled {new_name}"
                )
                utils.log("\tAll associated files will be referred to with this name")
                name = new_name
            # Save the data for this line.
            name_set.add(name)
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
    name_set: set[str] = set()
    for mol in suppl:
        # Convert mols to smiles. That's what the rest of the program is
        # designed to deal with.
        if mol:
            smiles = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
        else:
            utils.log(
                "\tWarning: Could not convert some SDF-formatted files to SMILES. Consider using an SMI (SMILES) file instead."
            )
            continue

        # A record with no atoms parses to a valid Mol, so it gets this far and
        # then has nothing to contribute. Drop it before the naming counters
        # advance, or it consumes an "untitled_N_molnum_N" slot and a name in
        # name_set that no output molecule ever uses, which shifts the names of
        # every record that follows.
        if smiles == "":
            utils.log("\tWarning: Skipping an SDF record that produced no atoms.")
            continue

        try:
            name = mol.GetProp("_Name")
        except Exception:
            name = ""

        # Handle unnamed ligands
        if not name:
            utils.log(
                f"\tUntitled ligand for the {mol_obj_counter} molecule in the input SDF"
            )
            name = f"untitled_{missing_name_counter}_molnum_{mol_obj_counter}"
            utils.log(f"\tNaming that ligand {name}")
            utils.log("\tAll associated files will be referred to with this name")
            missing_name_counter += 1

        # Handle duplicate ligands in same list.
        if name in name_set:
            utils.log(f"\nMultiple entries with the ligand name: {name}")
            # Keep incrementing the counter until the generated name is one
            # nothing else has claimed. An input that already contains, say,
            # "lig_copy_2" would otherwise collide with the name generated for
            # a duplicate "lig".
            duplicate_names[name] = duplicate_names.get(name, 1) + 1
            new_name = f"{name}_copy_{duplicate_names[name]}"
            while new_name in name_set:
                duplicate_names[name] += 1
                new_name = f"{name}_copy_{duplicate_names[name]}"
            name = new_name
            utils.log(
                f"\tThe version of the ligand for the {mol_obj_counter} molecule in the SDF file will be retitled {name}"
            )
            utils.log("\tAll associated files will be referred to with this name")
        mol_obj_counter += 1
        name_set.add(name)

        # SDF files may also contain properties. Get those as well.
        try:
            properties = mol.GetPropsAsDict()
        except Exception:
            properties = {}

        data.append((smiles, name, properties))

    return data
