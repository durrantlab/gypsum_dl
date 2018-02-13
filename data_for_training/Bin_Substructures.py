"""Bin_Substructures.py

This is a script for generating a substructure file for use in the 
AddHydrogens step of gypsum. This is an approximation for actual protonation at
best, and when a proper way of calculating protonation at various pHs is found
this should be replaced.
"""
import sys

from rdkit import Chem

import gypsum.Steps.SMILES.AddHydrogens as AddHs


def main(substructures, trainingData, threshold):
    """
    Main run script. Opens the given substructure folder, finds matching
    smiles, and calculates the pKa bin, and outputs it in the appropriate
    format for the AddHydrogens step.
    """

    # We first need to extract the substructures and their sites for protonation.
    # This has the inheirent assumption that explicit hydrogens in the SMART string
    # are the locations where we may protonate or deprotonate.
    # We do not have a way to mark those sites currently. We may wish to rethink
    # that at some later point.

    subs = []
    with open(substructures, 'r') as substruct:
        for line in substruct:
            line = line.strip()
            if line is not "":
                name, smart = line.split()
                smart_mol = Chem.MolFromSmarts(smart)
                prot = []
                for atom in smart_mol.GetAtoms():
                    if atom.GetAtomicNum() == 1:
                        for neighbor in atom.GetNeighbors():
                            prot.append(neighbor.GetIdx())
                subs.append((name, smart_mol, smart, prot))

    # We want to load the training data and assign its pKa to all matching
    # substructures
    pka_dict = {}

    with open(trainingData, 'r') as data:
        for line in data:
            line = line.strip()
            if line is not "":
                smile, pka = line.split('\t')
                mol = Chem.AddHs(Chem.MolFromSmiles(smile))
                matches = get_matches(mol, subs)

                for item in matches:
                    if item not in pka_dict.keys():
                        pka_dict[item] = [pka]
                    else:
                        pka_dict[item].append(pka)


    # Given our matched pKa values, we want to check if more than a certain
    # percentage, our threshold, are above 9 or below 5.
    bin_dict = {}

    for key in pka_dict:
        values = pka_dict[key]
        floats = [float(x) for x in values]

        total = float(len(floats))
        sub_5 = float(len([x for x in floats if x < 5.0]))
        sup_9 = float(len([x for x in floats if x > 9.0]))

        if (sub_5 / total) >= threshold:
            bin_dict[key] = "acid"
        elif (sup_9 / total) >= threshold:
            bin_dict[key] = "base"
        else:
            bin_dict[key] = "neutral"


    # Now we want to format and write our binned substructures to output
    output = []

    for item in subs:
        name, _, smart, prot = item
        if name in bin_dict.keys():
            category = bin_dict[name]
        else:
            category = "neutral" # We do not have enough information to bin the category

        temp = [name, smart]
        for site in prot:
            temp.append(str(site))
            temp.append(category)

        outline = "\t".join(temp)
        output.append(outline)

    with open("site_structures.smarts",'w') as out_file:
        out_file.write("\n".join(output))


def get_matches(mol, subs):
    """
    Given a molecule and a list of subgroups, return the subgroups that match
    where subgroups do not overlap.
    """
    AddHs.UnprotectMolecule(mol)
    mol_matches = []

    for item in subs:
        name, smart, smart_str, prot = item
        smart = Chem.MolFromSmarts(smart_str)
        #print(name, smart, smart_str)
        #print(mol.HasSubstructMatch(smart))
        if mol.HasSubstructMatch(smart):
            #print(name+"\t"+test_smile+"\t"+smart_str)
            matches = AddHs.UnprotectedMatches(mol, smart)
            for match in matches:
                mol_matches.append(name)
                AddHs.ProtectMolecule(mol, match)

    return mol_matches


if __name__ == "__main__":
    substructures = sys.argv[1] # "substructures-trimmed.smarts"
    trainingData = sys.argv[2]
    threshold = sys.argv[3]

    main(substructures, trainingData, threshold)
