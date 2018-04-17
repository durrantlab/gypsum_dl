import __future__

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)


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
    for line in open(filename):
        line = line.strip()
        if line != "":
            chunks = line.split()
            smiles = chunks[0]
            name = " ".join(chunks[1:])
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

        try:
            properties = mol.GetPropsAsDict()
        except:
            properties = {}

        if smiles != "":
            data.append((smiles, name, properties))
        
    return data
        

