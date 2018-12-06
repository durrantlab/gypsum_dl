import __future__

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")


def load_smiles_file(filename):
    """
    Loads a smiles file.

    :param str filename: The filename.

    :returns: A list of tuples, (SMILES, Name).
    :rtype: :class:`list[(str, str)]`
    """

    # The smiles file is the smiles string, separated by white space, followed
    # by the name.
    data = []        
    duplicate_names = {}
    missing_name_counter = 0 
    line_counter = 0
    name_list = []
    for line in open(filename):
        line = line.strip()
        if line != "":
            chunks = line.split()
            smiles = chunks[0]
            name = " ".join(chunks[1:])

            # Handle unnamed ligands
            if name == "":
                print("\tUntitled ligand on line {}".format(line_counter))
                name = "untitled_lig_{}_line_{}".format(missing_name_counter,line_counter)
                print("\tNaming that ligand {}".format(name))
                print("\tAll associated files will be refered to with this name")
                missing_name_counter += 1

            # Handle duplicate ligands in same list
            if name in name_list:
                # if multiple names 
                if name in list(duplicate_names.keys()):
                    duplicate_names[name] = duplicate_names[name] + 1

                    new_name = "{}_copy_{}".format(name, duplicate_names[name])
                    print("\nMultiple entries with the ligand name: {}".format(name))
                    print("\tThe veresion of the ligand on line {} will be retitled {}".format(line_counter,new_name))
                    print("\tAll associated files will be refered to with this name")
                    name = new_name
                else:
                    duplicate_names[name] = 2
                    new_name = "{}_copy_{}".format(name, duplicate_names[name])
                    print("\nMultiple entries with the ligand name: {}".format(name))
                    print("\tThe veresion of the ligand on line {} will be retitled {}".format(line_counter,new_name))
                    print("\tAll associated files will be refered to with this name")
                    name = new_name
            
            name_list.append(name)
            line_counter +=1
            data.append((smiles, name, {}))
    return data


def load_sdf_file(filename):
    """
    Loads an sdf file.

    :param str filename: The filename.

    :returns: A list of tuples, (SMILES, Name).
    :rtype: :class:`list[(str, str)]`
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
        smiles = Chem.MolToSmiles(
            mol, isomericSmiles=True, canonical=True
        )

        try:
            name = mol.GetProp("_Name")
        except:
            name = ""

        # Handle unnamed ligands
        if name == "":
            print("\tUntitled ligand for the {} molecule in the input SDF".format(mol_obj_counter))
            name = "untitled_lig_{}_molnum_{}".format(missing_name_counter,mol_obj_counter)
            print("\tNaming that ligand {}".format(name))
            print("\tAll associated files will be refered to with this name")
            missing_name_counter += 1

            # Handle duplicate ligands in same list
            if name in name_list:
                # if multiple names 
                if name in list(duplicate_names.keys()):
                    duplicate_names[name] = duplicate_names[name] + 1

                    new_name = "{}_copy_{}".format(name, duplicate_names[name])
                    print("\nMultiple entries with the ligand name: {}".format(name))
                    print("\tThe veresion of the ligand for the {} molecule in the SDF file will be retitled {}".format(mol_obj_counter,new_name))
                    print("\tAll associated files will be refered to with this name")
                    name = new_name
                else:
                    duplicate_names[name] = 2
                    new_name = "{}_copy_{}".format(name, duplicate_names[name])
                    print("\nMultiple entries with the ligand name: {}".format(name))
                    print("\tThe veresion of the ligand for the {} molecule in the SDF file will be retitled {}".format(mol_obj_counter,new_name))
                    print("\tAll associated files will be refered to with this name")
                    name = new_name

            mol_obj_counter += 1
            name_list.append(name)
        try:
            properties = mol.GetPropsAsDict()
        except:
            properties = {}

        if smiles != "":
            data.append((smiles, name, properties))
        
    return data
        

