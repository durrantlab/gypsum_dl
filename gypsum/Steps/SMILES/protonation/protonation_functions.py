"""
This script identifies and enumerates the possible protonation sites of SMILES strings.
"""
import copy
import os

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def protonate(args):
    """
    Protonates a set of molecules as given by the user inputs.
    """
    args = clean_args(args)
    subs = load_protonation_substructs(args["min_ph"], args["max_ph"], args["st_dev"])
    smiles = args["smiles"]
    data = args["data"]

    output = []
    for i, smi in enumerate(smiles):
        tag = " ".join(data[i])
        sites = get_protonation_sites(smi, subs)
        new_smis = [smi]
        for site in sites:
            new_smis = protonate_site(new_smis, site)
        new_lines = [x + '\t' + tag for x in new_smis]
        output.extend(new_lines)

    return output

def clean_args(args):
    """
    Cleans and normalizes input parameters
    """
    defaults = {'min_ph' : 6.4,
                'max_ph' : 8.4,
                'st_dev' : 1.0}

    for key in defaults:
        if key not in args:
            args[key] = defaults[key]

    keys = list(args.keys())
    for key in keys:
        if args[key] is None:
            del args[key]

    if "smiles" in args:
        if isinstance(args["smiles"], str):
            splits = args["smiles"].strip().split()
            args["smiles"] = [splits[0]]
            args["data"] = [splits[1:]]
    elif "smiles_file" in args:
        args["smiles"], args["data"] = load_files(args["smiles_file"])
    else:
        raise Exception("Error: No SMILES in params.")


    #print("START")
    mol_str_list = []
    for smiles_str in args["smiles"]:
        
        # Convert from SMILES string to RDKIT Mol
        # Filter if failed.
        mol = convert_smiles_str_to_mol(smiles_str)
        if mol is None:
            continue
                
        # Handle nuetralizing the molecules
        # Filter if failed.
        mol = neutralize_mol(mol)
        if mol is None:
            continue
        
        try:
            mol = Chem.RemoveHs(mol)
        except:
            continue

        if mol is None:
            continue

        new_mol_string = Chem.MolToSmiles(mol)
        mol_str_list.append(new_mol_string)


    args["smiles"] = [x for x in mol_str_list]

    return args

def neutralize_mol(mol):
    """
    All molecules need to be neuralized to the extent possible. The user should not be
    allowed to specify the valence of the atoms in most cases.
    """

    # Initialize some variables
    # To handle O- bonded to only one atom (add hydrogen).
    deprot_oxy = Chem.MolFromSmarts('[Ov1-1:1]')
    prot_oxy_rxn = None

    # To handle N+ bonded to a hydrogen (remove hydrogen).
    prot_nitrogen = Chem.MolFromSmarts('[#7v4+1:1]-[H]')
    deprot_nitrogen_rxn = None

    # To handle O- bonded to two atoms. Should not be Negative.
    wrong_prot_oxy = Chem.MolFromSmarts('[Ov2-:1]')
    wrong_prot_oxy_rxn = None

    # To handle N+ bonded to three atoms. Should not be positive.
    wrong_prot_nitrogen = Chem.MolFromSmarts('[#7v3+1:1]')
    wrong_prot_nitrogen_rxn = None

    # To handle N- Bonded to two atoms. Add hydrogen.
    wrong_prot_nitrogen2 = Chem.MolFromSmarts('[#7v2-1:1]')
    wrong_prot_nitrogen_rxn2 = None

    # Add hydrogens (respects valence, so incomplete).
    #Chem.calcImplicitValence(mol)
    mol.UpdatePropertyCache(strict=False)
    mol = Chem.AddHs(mol)

    while True:  # Keep going until all these issues have been resolved.
        rxn = None  # The reaction to perform.

        # Negative oxygen atom bonded to only one atom? Add hydrogen.
        if mol.HasSubstructMatch(deprot_oxy):
            if prot_oxy_rxn is None:
                prot_oxy_rxn = AllChem.ReactionFromSmarts('[Ov1-1:1]>>[Ov2+0:1]-[H]')
            rxn = prot_oxy_rxn

        # Positive, protonated nitrogen should be deprotonated
        elif mol.HasSubstructMatch(prot_nitrogen):
            if deprot_nitrogen_rxn is None:
                deprot_nitrogen_rxn = AllChem.ReactionFromSmarts('[#7v4+1:1]-[H]>>[#7v3+0:1]')
            rxn = deprot_nitrogen_rxn

        # Oxygen bonded to two atoms shouldn't be negative. I'm not so sure this could ever happen.
        elif mol.HasSubstructMatch(wrong_prot_oxy):
            if wrong_prot_oxy_rxn is None:
                wrong_prot_oxy_rxn = AllChem.ReactionFromSmarts('[Ov2-:1]>>[Ov2+0:1]')
            rxn = wrong_prot_oxy_rxn

        # Nitrogen bonded to three atoms shouldn't be pvesitie
        elif mol.HasSubstructMatch(wrong_prot_nitrogen):
            if wrong_prot_nitrogen_rxn is None:
                wrong_prot_nitrogen_rxn = AllChem.ReactionFromSmarts('[#7v3+1:1]>>[#7v3+0:1]')
            rxn = wrong_prot_nitrogen_rxn

        # Nitrogen bonded to two atoms shouldn't be negative. Need to add hydrogen.
        elif mol.HasSubstructMatch(wrong_prot_nitrogen2):
            if wrong_prot_nitrogen_rxn2 is None:
                wrong_prot_nitrogen_rxn2 = AllChem.ReactionFromSmarts('[#7v2-1:1]>>[#7+0:1]-[H]')
            rxn = wrong_prot_nitrogen_rxn2

        # Perform the reaction if necessary
        if rxn is None:  # No reaction left, so break out of while loop.
            break
        else:
            mol = rxn.RunReactants((mol,))[0][0]
            mol.UpdatePropertyCache(strict=False)  # Update valences

    # The mols have been altered from the reactions described above, we need to resanitize them.
    # Make sure aromatic rings are shown as such
    # This catches all RDKit Errors. without the catchError and sanitizeOps
    # the Chem.SanitizeMol can crash the program.
    sanitize_string =  Chem.SanitizeMol(mol, sanitizeOps = rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL, catchErrors = True)
    if sanitize_string.name == "SANITIZE_NONE": 
        
        return mol
    else:
        # Some Sanitation error occured.
        return None

def load_files(smile_file):
    """
    Loads smiles from file.
    """
    smiles = []
    data = []
    with open(smile_file, 'r') as smis:
        for line in smis:
            splits = line.split()
            if len(splits) != 0:
                smiles.append(splits[0])
                data.append(splits[1:])
    return smiles, data

def load_protonation_substructs(min_ph=6.4, max_ph=8.4, pka_std_range=1):
    """
    A pre-calculated list of R-groups with protonation sites, with their likely
    pKa bins.
    """
    subs = []
    pwd = os.path.dirname(__file__)
    site_structures_file = "{}/{}".format(pwd,"site_substructures.smarts")
    with open(site_structures_file, 'r') as substruct:
        for line in substruct:
            line = line.strip()
            sub = {}
            if line is not "":
                splits = line.split()
                sub["name"] = splits[0]
                sub["smart"] = splits[1]
                sub["mol"] = Chem.MolFromSmarts(sub["smart"])

                #NEED TO DIVIDE THIS BY 3s
                pka_ranges = [splits[i:i+3] for i in range(2, len(splits)-1, 3)]

                prot = []
                for pka_range in pka_ranges:
                    site = pka_range[0]
                    std = float(pka_range[2]) * pka_std_range
                    mean = float(pka_range[1])
                    protonation_state = define_protonation_state(mean, std, min_ph, \
                        max_ph)

                    prot.append([site, protonation_state])

                sub["prot"] = prot
                subs.append(sub)
    return subs

def define_protonation_state(mean, std, min_ph, max_ph):
    """
    Updates the substructure definitions to include the protonation state based on the user-given
    pH range. The size of the pKa range is also based on the number of standard deviations to be
    considered by the user param.
    """
    min_pka = mean - std
    max_pka = mean + std

    # This needs to be reassigned, and 'ERROR' should never make it past the next set of checks.
    protonation_state = 'ERROR'

    if min_pka <= max_ph and min_ph <= max_pka:
        protonation_state = 'BOTH'
    elif mean > max_ph:
        protonation_state = 'PROTONATED'
    elif mean < min_ph:
        protonation_state = 'DEPROTONATED'

    # We are error handling here
    if protonation_state == 'ERROR':
        raise Exception("HORRIBLE NONSENSE HAS OCCURED.")

    return protonation_state


###
# We need to identify and mark groups that have been matched with a substructure.
###

def unprotect_molecule(mol):
    """
    Sets the protected property on all atoms to 0. This also creates the property
    for new molecules.
    """
    for atom in mol.GetAtoms():
        atom.SetProp('_protected', '0')

def protect_molecule(mol, match):
    """
    Given a 'match', a list of molecules idx's, we set the protected status of each
    atom to 1. This will prevent any matches using that atom in the future.
    """
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        atom.SetProp('_protected', '1')

def get_unprotected_matches(mol, substruct):
    """
    Finds substructure matches with atoms that have not been protected.
    Returns list of matches, each match a list of atom idxs.
    """
    matches = mol.GetSubstructMatches(substruct)
    unprotected_matches = []
    for match in matches:
        if is_match_unprotected(mol, match):
            unprotected_matches.append(match)
    return unprotected_matches

def is_match_unprotected(mol, match):
    """
    Checks a molecule to see if the substructure match
    contains any protected atoms.
    """
    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        protected = atom.GetProp("_protected")
        if protected == "1":
            return False
    return True

def neutralize_molecule(mol):
    """
    Neutralize things. Maybe?
    """
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)

def get_protonation_sites(smi, subs):
    """
    For a single molecule, find all possible matches in the protonation R-group list,
    subs. Items that are higher on the list will be matched first, to the exclusion of
    later items.
    Returns a list of protonation sites and their pKa bin. ('Acid', 'Neutral', or 'Base')
    """
    # Convert the Smiles string (smi) to an RDKit Mol Obj
    mol = convert_smiles_str_to_mol(smi)

    # Check Conversion worked
    if mol is None:
        print("ERROR:   ",smi)
        return []
        
    # Try to Add hydrogens. if failed return []
    try:
        mol =  Chem.AddHs(mol)
    except:
        print("ERROR:   ",smi)
        return []

    # Check adding Hs worked
    if mol is None:
        print("ERROR:   ",smi)
        return []

    unprotect_molecule(mol)
    protonation_sites = []

    for item in subs:
        smart = item['mol']
        if mol.HasSubstructMatch(smart):
            matches = get_unprotected_matches(mol, smart)
            prot = item['prot']
            for match in matches:
                # We want to move the site from being relative to the
                # substructure, to the index on the main molecule.
                for site in prot:
                    proton = int(site[0])
                    category = site[1]
                    new_site = (match[proton], category)
                    protonation_sites.append(new_site)
                protect_molecule(mol, match)
    return protonation_sites

def protonate_site(smis, site):
    """
    Given a list of smis, we protonate the site.
    """
    # Decouple the atom index and its charge from the site tuple
    idx, charge = site

    # Initialize the output list
    output_smis = []

    charge_dict = {"DEPROTONATED": [-1],
                   "PROTONATED": [0],
                   "BOTH": [-1, 0]}

    charges = charge_dict[charge]

    output_smis = set_protonation_charge(smis, idx, charges)

    return output_smis

def set_protonation_charge(smis, idx, charges):
    """
    Sets the atomic charge on a particular site for a set of SMILES.
    """
    # Sets up the output list and the Nitrogen charge
    output = []

    for charge in charges:
        # The charge for Nitrogens is 1 higher than others
        nitro_charge = charge + 1

        for smi in smis:
            
            # Convert smilesstring (smi) into a RDKit Mol
            mol = convert_smiles_str_to_mol(smi)
             
            # Check that the conversion worked, skip if it fails
            if mol is None:
                continue    
                            
            atom = mol.GetAtomWithIdx(idx)

            # Assign the protonation charge, with special care for Nitrogens
            element = atom.GetAtomicNum()
            if element == 7:
                atom.SetFormalCharge(nitro_charge)
            else:
                atom.SetFormalCharge(charge)

            # Convert back to SMILE and add to output
            out_smile = Chem.MolToSmiles(mol, isomericSmiles=True,canonical=True)
            output.append(out_smile)

    return output

def convert_smiles_str_to_mol(smiles_str):
    """
    Given a SMILES string check that it is actually a string and not a None
    then try to convert it to an RDKit Mol Object.
    Return None if it is the wrong type or if it fails to convert to a Mol Obj
    Return the mol object if it converts.
    """
   
    if smiles_str is None or type(smiles_str) is not str:
        return None

    # Check that there are no type errors, ie Nones or non-string
    # A non-string type will cause RDKit to hard crash
    try:
        mol = Chem.MolFromSmiles(smiles_str)
    except:
        return None

    # Check that there are None type errors
    # Chem.MolFromSmiles has sanitize on which means if there is
    #   even a small error in the SMILES (kekulize, nitrogen charge...) then mol=None. 
    #       ie. Chem.MolFromSmiles("C[N]=[N]=[N]") = None
    #           this is an example of an nitrogen charge error.
    #   It is cased in a try statement to be overly cautious.
     
    if mol is None:

        return None

    else:
        return mol

