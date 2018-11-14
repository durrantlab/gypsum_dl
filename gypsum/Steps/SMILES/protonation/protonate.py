# Copyright 2018 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This script identifies and enumerates the possible protonation sites of SMILES
strings.
"""

import copy
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import argparse

def main():
    """The main definition run when you call the script from the commandline."""

    parser = get_args()
    args = vars(parser.parse_args())

    if args["test"]:
        # Run tests.
        test()
    else:
        # Run protonation
        output = protonate(args)

        if 'output_file' in args:
            with open(args['output_file'], 'w') as file:
                file.write("\n".join(output))
        else:
            for out in output:
                print(out)

def get_args():
    """Gets the arguments from the command line.

    :return: A parser object.
    """

    parser = argparse.ArgumentParser(description="Protonates small moleucles.")
    parser.add_argument('--min_ph', metavar='MIN', type=float, default=6.4,
                        help='minimum pH to consider')
    parser.add_argument('--max_ph', metavar='MAX', type=float, default=8.4,
                        help='maximum pH to consider')
    parser.add_argument('--pka_precision', metavar='PRE', type=float, default=1.0,
                        help='pKa precision factor (number of standard devations)')
    parser.add_argument('--smiles', type=str,
                        help='SMILES string to protonate')
    parser.add_argument('--smiles_file', type=str,
                        help='file that contains SMILES strings to protonate')
    parser.add_argument('--output_file', type=str,
                        help='output file to write protonated SMILES (optional)')
    parser.add_argument('--label_states', action="store_true",
                        help='label protonated SMILES with target state ' + \
                             '("DEPROTONATED", "PROTONATED", or "BOTH").')
    parser.add_argument('--test', action="store_true",
                        help='run unit tests (for debugging)')

    return parser

def protonate(args):
    """Protonates a set of molecules as given by the user inputs.

    :param dict args: A dictionary containing the arguments.
    :return: A list of the protonated smiles strings.
    """

    args = clean_args(args)
    subs = load_protonation_substructs_calc_state_for_ph(
        args["min_ph"], args["max_ph"], args["pka_precision"]
    )
    smiles = args["smiles"]

    # Note that args["data"] includes everything on SMILES line but the SMILES
    # string itself (e.g., the molecule name). It is set in clean_args().
    data = args["data"]

    output = []
    for i, smi in enumerate(smiles):
        if smi.startswith("NONE|"):
            print("ERROR: Skipping poorly formed SMILES string: " + smi[5:] + "\t" + " ".join(data[i]))
            continue

        # Collect the data associated with this smiles (e.g., the molecule
        # name).
        tag = " ".join(data[i])

        # sites is a list of (atom index, "PROTONATED|DEPROTONATED|BOTH").
        # Note that the second entry indicates what state the site SHOULD be
        # in (not the one it IS in per the SMILES string). It's calculated
        # based on the probablistic distributions obtained during training.
        sites = get_prot_sites_and_target_states(smi, subs)

        new_smis = [smi]
        for site in sites:
            # Make a new smiles with the correct protonation state. Note that
            # new_smis is a growing list. This is how multiple protonation
            # sites are handled.

            # new_smis_to_perhaps_add = protonate_site(new_smis, site)
            new_smis = protonate_site(new_smis, site)

            # Only add new smiles if not already in the list.
            # for s in new_smis_to_perhaps_add:
                # if not s in new_smis:
                    # new_smis.append(s)

        # In some cases, the script might generate redundant molecules.
        # Phosphonates, when the pH is between the two pKa values and the
        # stdev value is big enough, for example, will generate two identical
        # BOTH states. Let's remove this redundancy.
        new_smis = list(set(new_smis))

        # If the user wants to see the target states, add those
        # to the ends of each line.
        if args["label_states"]:
            states = '\t'.join([x[1] for x in sites])
            new_lines = [x + "\t" + tag + "\t" + states for x in new_smis]
        else:
            new_lines = [x + "\t" + tag for x in new_smis]

        output.extend(new_lines)

    return output

def clean_args(args):
    """Cleans and normalizes input parameters

    :param dict args: A dictionary containing the arguments.
    :raises Exception: No SMILES in params.
    :return: A dict of the arguments, now fixed.
    """

    defaults = {'min_ph' : 6.4,
                'max_ph' : 8.4,
                'pka_precision' : 1.5,
                'label_states' : False,
                'test' : False}

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

    mol_str_list = []
    for i, smiles_str in enumerate(args["smiles"]):

        # Convert from SMILES string to RDKIT Mol
        # Filter if failed.

        mol = convert_smiles_str_to_mol(smiles_str)
        if mol is None:
            mol_str_list.append("NONE|" + smiles_str)
            continue

        # Handle nuetralizing the molecules
        # Filter if failed.
        mol = neutralize_mol(mol)
        if mol is None:
            mol_str_list.append("NONE|" + smiles_str)
            continue

        try:
            mol = Chem.RemoveHs(mol)
        except:
            mol_str_list.append("NONE|" + smiles_str)
            continue

        if mol is None:
            mol_str_list.append("NONE|" + smiles_str)
            continue

        new_mol_string = Chem.MolToSmiles(mol)
        mol_str_list.append(new_mol_string)

    args["smiles"] = [x for x in mol_str_list]

    return args

def neutralize_mol(mol):
    """All molecules need to be neuralized to the extent possible. The user
    should not be allowed to specify the valence of the atoms in most cases.

    :param rdkit.Chem.rdchem.Mol mol: The rdkit Mol objet to be neutralized.
    :return: The neutralized Mol object.
    """

    # Get the reaction data
    rxn_data = [
        ['[Ov1-1:1]', '[Ov2+0:1]-[H]'],  # To handle O- bonded to only one atom (add hydrogen).
        ['[#7v4+1:1]-[H]', '[#7v3+0:1]'],  # To handle N+ bonded to a hydrogen (remove hydrogen).
        ['[Ov2-:1]', '[Ov2+0:1]'],  # To handle O- bonded to two atoms. Should not be Negative.
        ['[#7v3+1:1]', '[#7v3+0:1]'],  # To handle N+ bonded to three atoms. Should not be positive.
        ['[#7v2-1:1]', '[#7+0:1]-[H]'],  # To handle N- Bonded to two atoms. Add hydrogen.
        # ['[N:1]=[N+0:2]=[N:3]-[H]', '[N:1]=[N+1:2]=[N+0:3]-[H]'],  # To handle bad azide. Must be protonated. (Now handled elsewhere, before SMILES converted to Mol object.)
        ['[H]-[N:1]-[N:2]#[N:3]', '[N:1]=[N+1:2]=[N:3]-[H]']  # To handle bad azide. R-N-N#N should be R-N=[N+]=N
    ]

    # Add substructures and reactions (initially none)
    for i, rxn_datum in enumerate(rxn_data):
        rxn_data[i].append(Chem.MolFromSmarts(rxn_datum[0]))
        rxn_data[i].append(None)

    # Add hydrogens (respects valence, so incomplete).
    # Chem.calcImplicitValence(mol)
    mol.UpdatePropertyCache(strict=False)
    mol = Chem.AddHs(mol)

    while True:  # Keep going until all these issues have been resolved.
        current_rxn = None  # The reaction to perform.
        current_rxn_str = None

        for i, rxn_datum in enumerate(rxn_data):
            reactant_smarts, product_smarts, substruct_match_mol, rxn_placeholder = rxn_datum
            if mol.HasSubstructMatch(substruct_match_mol):
                if rxn_placeholder is None:
                    current_rxn_str = reactant_smarts + '>>' + product_smarts
                    current_rxn = AllChem.ReactionFromSmarts(current_rxn_str)
                    rxn_data[i][3] = current_rxn  # Update the placeholder.
                else:
                    current_rxn = rxn_data[i][3]
                break

        # Perform the reaction if necessary
        if current_rxn is None:  # No reaction left, so break out of while loop.
            break
        else:
            mol = current_rxn.RunReactants((mol,))[0][0]
            mol.UpdatePropertyCache(strict=False)  # Update valences

    # The mols have been altered from the reactions described above, we need
    # to resanitize them. Make sure aromatic rings are shown as such This
    # catches all RDKit Errors. without the catchError and sanitizeOps the
    # Chem.SanitizeMol can crash the program.
    sanitize_string =  Chem.SanitizeMol(
        mol,
        sanitizeOps=rdkit.Chem.rdmolops.SanitizeFlags.SANITIZE_ALL,
        catchErrors = True
    )

    return mol if sanitize_string.name == "SANITIZE_NONE" else None

def load_files(smile_file):
    """Loads smiles from file.

    :param string smile_file: The smiles file name.
    :return: (list, list), the smiles strings and associated data,
             respectively.
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

def load_protonation_substructs_calc_state_for_ph(min_ph=6.4, max_ph=8.4, pka_std_range=1):
    """A pre-calculated list of R-groups with protonation sites, with their
    likely pKa bins.

    :param float min_ph:  The lower bound on the pH range, defaults to 6.4.
    :param float max_ph:  The upper bound on the pH range, defaults to 8.4.
    :param pka_std_range: Basically the precision (stdev from predicted pKa to
                          consider), defaults to 1.
    :return: A dict of the protonation substructions for the specified pH
            range.
    """

    subs = []
    pwd = os.path.dirname(os.path.realpath(__file__))

    site_structures_file = "{}/{}".format(pwd, "site_substructures.smarts")
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
                    protonation_state = define_protonation_state(
                        mean, std, min_ph, max_ph
                    )

                    prot.append([site, protonation_state])

                sub["prot_states_for_pH"] = prot
                subs.append(sub)
    return subs

def define_protonation_state(mean, std, min_ph, max_ph):
    """Updates the substructure definitions to include the protonation state
    based on the user-given pH range. The size of the pKa range is also based
    on the number of standard deviations to be considered by the user param.

    :param float mean:   The mean pKa.
    :param float std:    The precision (stdev).
    :param float min_ph: The min pH of the range.
    :param float max_ph: The max pH of the range.
    :raises Exception:   HORRIBLE NONSENSE HAS OCCURED.
    :return: A string describing the protonation state.
    """

    min_pka = mean - std
    max_pka = mean + std

    # This needs to be reassigned, and 'ERROR' should never make it past the
    # next set of checks.
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
# We need to identify and mark groups that have been matched with a
# substructure.
###

def unprotect_molecule(mol):
    """Sets the protected property on all atoms to 0. This also creates the
    property for new molecules.

    :param rdkit.Chem.rdchem.Mol mol: The rdkit Mol object.
    :type mol: The rdkit Mol object with atoms unprotected.
    """

    for atom in mol.GetAtoms():
        atom.SetProp('_protected', '0')

def protect_molecule(mol, match):
    """Given a 'match', a list of molecules idx's, we set the protected status
    of each atom to 1. This will prevent any matches using that atom in the
    future.

    :param rdkit.Chem.rdchem.Mol mol: The rdkit Mol object to protect.
    :param list match: A list of molecule idx's.
    """

    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        atom.SetProp('_protected', '1')

def get_unprotected_matches(mol, substruct):
    """Finds substructure matches with atoms that have not been protected.
    Returns list of matches, each match a list of atom idxs.

    :param rdkit.Chem.rdchem.Mol mol: The Mol object to consider.
    :param string substruct: The SMARTS string of the substructure ot match.
    :return: A list of the matches. Each match is itself a list of atom idxs.
    """

    matches = mol.GetSubstructMatches(substruct)
    unprotected_matches = []
    for match in matches:
        if is_match_unprotected(mol, match):
            unprotected_matches.append(match)
    return unprotected_matches

def is_match_unprotected(mol, match):
    """Checks a molecule to see if the substructure match contains any
    protected atoms.

    :param rdkit.Chem.rdchem.Mol mol: The Mol object to check.
    :param list match: The match to check.
    :return: A boolean, whether the match is present or not.
    """

    for idx in match:
        atom = mol.GetAtomWithIdx(idx)
        protected = atom.GetProp("_protected")
        if protected == "1":
            return False
    return True

def neutralize_molecule(mol):
    """Neutralize things. Maybe?

    :param rdkit.Chem.rdchem.Mol mol: The Mol object to conside.r
    """

    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)

def get_prot_sites_and_target_states(smi, subs):
    """For a single molecule, find all possible matches in the protonation
    R-group list, subs. Items that are higher on the list will be matched
    first, to the exclusion of later items.

    :param string smi: A SMILES string.
    :param list subs: Substructure information.
    :return: A list of protonation sites and their pKa bin. ('PROTONATED',
        'BOTH', or  'DEPROTONATED')
    """

    # Convert the Smiles string (smi) to an RDKit Mol Obj
    mol = convert_smiles_str_to_mol(smi)

    # Check Conversion worked
    if mol is None:
        print("ERROR:   ", smi)
        return []

    # Try to Add hydrogens. if failed return []
    try:
        mol =  Chem.AddHs(mol)
    except:
        print("ERROR:   ", smi)
        return []

    # Check adding Hs worked
    if mol is None:
        print("ERROR:   ", smi)
        return []

    unprotect_molecule(mol)
    protonation_sites = []

    for item in subs:
        smart = item["mol"]
        if mol.HasSubstructMatch(smart):
            matches = get_unprotected_matches(mol, smart)
            prot = item["prot_states_for_pH"]
            for match in matches:
                # We want to move the site from being relative to the
                # substructure, to the index on the main molecule.
                for site in prot:
                    proton = int(site[0])
                    category = site[1]
                    new_site = (match[proton], category, item["name"])

                    if not new_site in protonation_sites:
                        # Because sites must be unique.
                        protonation_sites.append(new_site)

                protect_molecule(mol, match)

    return protonation_sites

def protonate_site(smis, site):
    """Given a list of SMILES strings, we protonate the site.

    :param list smis:  The list of SMILES strings.
    :param tuple site: Information about the protonation site.
                       (idx, target_prot_state, prot_site_name)
    :return: A list of the appropriately protonated SMILES.
    """

    # Decouple the atom index and its target protonation state from the site
    # tuple
    idx, target_prot_state, prot_site_name = site

    # Initialize the output list
    output_smis = []

    state_to_charge = {"DEPROTONATED": [-1],
                       "PROTONATED": [0],
                       "BOTH": [-1, 0]}

    charges = state_to_charge[target_prot_state]

    # Now make the actual smiles match the target protonation state.
    output_smis = set_protonation_charge(smis, idx, charges, prot_site_name)

    return output_smis

def set_protonation_charge(smis, idx, charges, prot_site_name):
    """Sets the atomic charge on a particular site for a set of SMILES.

    :param list smis:             A list of the SMILES strings.
    :param int idx:               The index of the atom to consider.
    :param list charges:          A list of the charges (ints) to assign at
                                  this site.
    :param string prot_site_name: The name of the protonation site.
    :return: A list of the processed SMILES strings.
    """

    # Sets up the output list and the Nitrogen charge
    output = []

    for charge in charges:
        # The charge for Nitrogens is 1 higher than others (i.e., protonated
        # state is positively charged).
        nitro_charge = charge + 1

        # But there are a few nitrogen moieties where the acidic group is the
        # neutral one. Amides are a good example. I gave some thought re. how
        # to best flag these. I decided that those nitrogen-containing
        # moieties where the acidic group is neutral (rather than positively
        # charged) will have "*" in the name.
        if "*" in prot_site_name:
            nitro_charge = nitro_charge - 1  # Undo what was done previously.

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
    """Given a SMILES string, check that it is actually a string and not a
    None. Then try to convert it to an RDKit Mol Object.

    :param string smiles_str: The SMILES string.
    :return: A rdkit.Chem.rdchem.Mol object, or None if it is the wrong type or
        if it fails to convert to a Mol Obj
    """

    if smiles_str is None or type(smiles_str) is not str:
        return None

    # Check that there are no type errors, ie Nones or non-string
    # A non-string type will cause RDKit to hard crash
    try:
        # Try to fix azides here. They are just tricky to deal with.
        smiles_str = smiles_str.replace("N=N=N", "N=[N+]=N")
        smiles_str = smiles_str.replace("NN#N", "N=[N+]=N")

        mol = Chem.MolFromSmiles(smiles_str)
    except:
        return None

    # Check that there are None type errors Chem.MolFromSmiles has sanitize on
    # which means if there is even a small error in the SMILES (kekulize,
    # nitrogen charge...) then mol=None. ie.
    # Chem.MolFromSmiles("C[N]=[N]=[N]") = None this is an example of an
    # nitrogen charge error. It is cased in a try statement to be overly
    # cautious.

    return None if mol is None else mol

def test():
    """Tests all the 38 groups."""

    smis = [
        # [input smiles, pka, protonated, deprotonated, category]
        ["C#CCO",                  "C#CCO",                     "C#CC[O-]",                 "Alcohol"],
        ["C(=O)N",                 "NC=O",                      "[NH-]C=O",                 "Amide"],
        ["CC(=O)NOC(C)=O",         "CC(=O)NOC(C)=O",            "CC(=O)[N-]OC(C)=O",        "Amide_electronegative"],
        ["COC(=N)N",               "COC(N)=[NH2+]",             "COC(=N)N",                 "AmidineGuanidine2"],
        ["Brc1ccc(C2NCCS2)cc1",    "Brc1ccc(C2[NH2+]CCS2)cc1",  "Brc1ccc(C2NCCS2)cc1",      "Amines_primary_secondary_tertiary"],
        ["CC(=O)[n+]1ccc(N)cc1",   "CC(=O)[n+]1ccc([NH3+])cc1", "CC(=O)[n+]1ccc(N)cc1",     "Anilines_primary"],
        ["CCNc1ccccc1",            "CC[NH2+]c1ccccc1",          "CCNc1ccccc1",              "Anilines_secondary"],
        ["Cc1ccccc1N(C)C",         "Cc1ccccc1[NH+](C)C",        "Cc1ccccc1N(C)C",           "Anilines_tertiary"],
        ["BrC1=CC2=C(C=C1)NC=C2",  "Brc1ccc2[nH+]ccc2c1",       "Brc1ccc2[nH]ccc2c1",       "Aromatic_nitrogen_protonated"],
        ["C-N=[N+]=[N@H]",         "CN=[N+]=N",                 "CN=[N+]=[N-]",             "Azide"],
        ["BrC(C(O)=O)CBr",         "O=C(O)C(Br)CBr",            "O=C([O-])C(Br)CBr",        "Carboxyl"],
        ["NC(NN=O)=N",             "NC(=[NH2+])NN=O",           "N=C(N)NN=O",               "AmidineGuanidine1"],
        ["C(F)(F)(F)C(=O)NC(=O)C", "CC(=O)NC(=O)C(F)(F)F",      "CC(=O)[N-]C(=O)C(F)(F)F",  "Imide"],
        ["O=C(C)NC(C)=O",          "CC(=O)NC(C)=O",             "CC(=O)[N-]C(C)=O",         "Imide2"],
        ["CC(C)(C)C(N(C)O)=O",     "CN(O)C(=O)C(C)(C)C",        "CN([O-])C(=O)C(C)(C)C",    "N-hydroxyamide"],
        ["C[N+](O)=O",             "C[N+](=O)O",                "C[N+](=O)[O-]",            "Nitro"],
        ["O=C1C=C(O)CC1",          "O=C1C=C(O)CC1",             "O=C1C=C([O-])CC1",         "O=C-C=C-OH"],
        ["C1CC1OO",                "OOC1CC1",                   "[O-]OC1CC1",               "Peroxide2"],
        ["C(=O)OO",                "O=COO",                     "O=CO[O-]",                 "Peroxide1"],
        ["Brc1cc(O)cc(Br)c1",      "Oc1cc(Br)cc(Br)c1",         "[O-]c1cc(Br)cc(Br)c1",     "Phenol"],
        ["CC(=O)c1ccc(S)cc1",      "CC(=O)c1ccc(S)cc1",         "CC(=O)c1ccc([S-])cc1",     "Phenyl_Thiol"],
        ["C=CCOc1ccc(C(=O)O)cc1",  "C=CCOc1ccc(C(=O)O)cc1",     "C=CCOc1ccc(C(=O)[O-])cc1", "Phenyl_carboxyl"],
        ["COP(=O)(O)OC",           "COP(=O)(O)OC",              "COP(=O)([O-])OC",          "Phosphate_diester"],
        ["CP(C)(=O)O",             "CP(C)(=O)O",                "CP(C)(=O)[O-]",            "Phosphinic_acid"],
        ["CC(C)OP(C)(=O)O",        "CC(C)OP(C)(=O)O",           "CC(C)OP(C)(=O)[O-]",       "Phosphonate_ester"],
        ["CC1(C)OC(=O)NC1=O",      "CC1(C)OC(=O)NC1=O",         "CC1(C)OC(=O)[N-]C1=O",     "Ringed_imide1"],
        ["O=C(N1)C=CC1=O",         "O=C1C=CC(=O)N1",            "O=C1C=CC(=O)[N-]1",        "Ringed_imide2"],
        ["O=S(OC)(O)=O",           "COS(=O)(=O)O",              "COS(=O)(=O)[O-]",          "Sulfate"],
        ["COc1ccc(S(=O)O)cc1",     "COc1ccc(S(=O)O)cc1",        "COc1ccc(S(=O)[O-])cc1",    "Sulfinic_acid"],
        ["CS(N)(=O)=O",            "CS(N)(=O)=O",               "CS([NH-])(=O)=O",          "Sulfonamide"],
        ["CC(=O)CSCCS(O)(=O)=O",   "CC(=O)CSCCS(=O)(=O)O",      "CC(=O)CSCCS(=O)(=O)[O-]",  "Sulfonate"],
        ["CC(=O)S",                "CC(=O)S",                   "CC(=O)[S-]",               "Thioic_acid"],
        ["C(C)(C)(C)(S)",          "CC(C)(C)S",                 "CC(C)(C)[S-]",             "Thiol"],
        ["Brc1cc[nH+]cc1",         "Brc1cc[nH+]cc1",            "Brc1ccncc1",               "Aromatic_nitrogen_unprotonated"],
        ["C=C(O)c1c(C)cc(C)cc1C",  "C=C(O)c1c(C)cc(C)cc1C",     "C=C([O-])c1c(C)cc(C)cc1C", "Vinyl_alcohol"],
        ["CC(=O)ON",               "CC(=O)O[NH3+]",             "CC(=O)ON",                 "Primary_hydroxyl_amine"]
    ]

    smis_phos = [
        ["O=P(O)(O)OCCCC", "CCCCOP(=O)(O)O", "CCCCOP(=O)([O-])O", "CCCCOP(=O)([O-])[O-]", "Phosphate"],
        ["CC(P(O)(O)=O)C", "CC(C)P(=O)(O)O", "CC(C)P(=O)([O-])O", "CC(C)P(=O)([O-])[O-]", "Phosphonate"]
    ]

    # Load the average pKa values.
    average_pkas = {l.split()[0].replace("*", ""):float(l.split()[3]) for l in open("site_substructures.smarts") if l.split()[0] not in ["Phosphate", "Phosphonate"]}
    average_pkas_phos = {l.split()[0].replace("*", ""):[float(l.split()[3]), float(l.split()[6])] for l in open("site_substructures.smarts") if l.split()[0] in ["Phosphate", "Phosphonate"]}

    print("")
    print("Running Tests")
    print("=============")
    print("")

    print("Very Acidic (pH -10000000)")
    print("--------------------------")
    print("")

    args = {
        "min_ph": -10000000,
        "max_ph": -10000000,
        "pka_precision": 0.5,
        "smiles": "",
        "label_states": True
    }

    for smi, protonated, deprotonated, category in smis:
        args["smiles"] = smi
        _test_check(args, [protonated], ["PROTONATED"])

    for smi, protonated, mix, deprotonated, category in smis_phos:
        args["smiles"] = smi
        _test_check(args, [protonated], ["PROTONATED"])

    args["min_ph"] = 10000000
    args["max_ph"] = 10000000

    print("")
    print("Very Basic (pH 10000000)")
    print("------------------------")
    print("")

    for smi, protonated, deprotonated, category in smis:
        args["smiles"] = smi
        _test_check(args, [deprotonated], ["DEPROTONATED"])

    for smi, protonated, mix, deprotonated, category in smis_phos:
        args["smiles"] = smi
        _test_check(args, [deprotonated], ["DEPROTONATED"])

    print("")
    print("pH is Category pKa")
    print("------------------")
    print("")

    for smi, protonated, deprotonated, category in smis:
        avg_pka = average_pkas[category]

        args["smiles"] = smi
        args["min_ph"] = avg_pka
        args["max_ph"] = avg_pka

        _test_check(args, [protonated, deprotonated], ["BOTH"])

    for smi, protonated, mix, deprotonated, category in smis_phos:
        args["smiles"] = smi

        avg_pka = average_pkas_phos[category][0]
        args["min_ph"] = avg_pka
        args["max_ph"] = avg_pka

        _test_check(args, [mix, protonated], ["BOTH"])

        avg_pka = average_pkas_phos[category][1]
        args["min_ph"] = avg_pka
        args["max_ph"] = avg_pka

        _test_check(args, [mix, deprotonated], ["DEPROTONATED", "DEPROTONATED"])

        avg_pka = 0.5 * (average_pkas_phos[category][0] + average_pkas_phos[category][1])
        args["min_ph"] = avg_pka
        args["max_ph"] = avg_pka
        args["pka_precision"] = 5  # Should give all three

        _test_check(args, [mix, deprotonated, protonated], ["BOTH", "BOTH"])

def _test_check(args, expected_output, labels):
    """Tests most ionizable groups. The ones that can only loose or gain a single proton.

    :param args: The arguments to pass to protonate()
    :param expected_output: A list of the expected SMILES-strings output.
    :param labels: The labels. A list containing combo of BOTH, PROTONATED,
                   DEPROTONATED.
    :raises Exception: Wrong number of states produced.
    :raises Exception: Unexpected output SMILES.
    :raises Exception: Wrong labels.
    """

    output = protonate(args)
    output = [o.split() for o in output]

    num_states = len(expected_output)

    if (len(output) != num_states):
        raise Exception(
            args["smiles"][0] + " should have " + str(num_states) + " states at at pH " + str(args["min_ph"]) + ": " + \
            str(output)
        )

    if (len(set([l[0] for l in output]) - set(expected_output)) != 0):
        raise Exception(
            args["smiles"][0] + " is not " + " AND ".join(expected_output) + " at pH " + str(args["min_ph"]) + " - " + str(args["max_ph"]) + "; it is " + " AND ".join([l[0] for l in output])
        )

    if (len(set([l[1] for l in output]) - set(labels)) != 0):
        raise Exception(
            args["smiles"][0] + " not labeled as " + " AND ".join(labels) + "; it is " + " AND ".join([l[1] for l in output])
        )

    ph_range = sorted(list(set([args["min_ph"], args["max_ph"]])))
    ph_range_str = "(" + " - ".join("{0:.2f}".format(n) for n in ph_range) + ")"
    print("(CORRECT) " + ph_range_str.ljust(10) + " " + args["smiles"][0] + " => " + " AND ".join([l[0] for l in output]))

if __name__ == "__main__":
    main()
