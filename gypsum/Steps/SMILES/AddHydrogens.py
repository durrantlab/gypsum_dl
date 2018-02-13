"""
This module is made to identify and enumerate the possible protonation sites of molecules.
"""

# TODO: Need to have categories, potentially for each site in an R-group?

import copy

from rdkit import Chem
import gypsum.Multiprocess as mp
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol

import os

def load_protonation_substructs():
    """
    A pre-calculated list of R-groups with protonation sites, with their likely pKa bins.
    """
    subs = []
    file = os.path.join(os.path.dirname(__file__),"site_structures.smarts")
    with open(file, 'r') as substruct:
        for line in substruct:
            line = line.strip()
            sub = {}
            if line is not "":
                splits = line.split()
                name = splits[0]
                smart = splits[1]

                sub["name"] = name
                sub["smart"] = smart

                mol = Chem.MolFromSmarts(smart)
                sub["mol"] = mol

                # This is going to split the remaining
                prot = [[int(splits[i]), splits[i+1]] for i in range(2, len(splits)-1, 2)]

                #!!! This takes any explicit hydrogen as a site for protonation.
                # We may find that we need to revisit this, as Nitrogen causes issues
                # with RDKit's kekulization when it has 4 bonds in aromatic rings.
                sub["prot"] = prot
                subs.append(sub)
    return subs


###
# We need to identify and mark groups that have been matched with a substructure.
###

def remove_explicit_Hs(mol):
    for atom in mol.GetAtoms():
        atom.SetNumExplicitHs(0)

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
        keep_flag = True
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            protected = atom.GetProp("_protected")
            if protected == "1":
                keep_flag = False
        if keep_flag:
            unprotected_matches.append(match)
            protect_molecule(mol, match)
    return unprotected_matches

def get_protonation_sites(mol, subs):
    """
    For a single molecule, find all possible matches in the protonation R-group list,
    subs. Items that are higher on the list will be matched first, to the exclusion of
    later items.
    Returns a list of protonation sites and their pKa bin. ('Acid', 'Neutral', or 'Base')
    """
    unprotect_molecule(mol)
    remove_explicit_Hs(mol)
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
                    proton = site[0]
                    category = site[1]
                    new_site = (match[proton], category)
                    protonation_sites.append(new_site)
                protect_molecule(mol, match)

    return protonation_sites

def protonate_site(mols, site):
    """
    Protonates (or deprotonates) a list of smis at a given site.
    """
    # We get the index to modify and the charges to apply
    idx, charge = site

    # Create lists of things to positively or negatively charge.
    # Neutral gets added to both lists
    protonate = []
    deprotonate = []

    if charge == "acid" or charge == "neutral":
        deprotonate = copy.deepcopy(mols)
    if charge == "base" or charge == "neutral":
        protonate = copy.deepcopy(mols)

    output_mols = []

    for mol in deprotonate:
        mol = Chem.RemoveHs(mol)
        atom = mol.GetAtomWithIdx(idx)
        atom.SetNumExplicitHs(0)
        element = atom.GetAtomicNum()
        if element == 7:
            atom.SetFormalCharge(0)
        else:
            #print("Minus")
            atom.SetFormalCharge(-1)
        output_mols.append(mol)

    for mol in protonate:
        mol = Chem.RemoveHs(mol)
        atom = mol.GetAtomWithIdx(idx)
        element = atom.GetAtomicNum()
        if element == 7:
            atom.SetFormalCharge(+1)
        else:
            #print("Neutral")
            atom.SetFormalCharge(0)
        output_mols.append(mol)

    return output_mols

def add_hydrogens(self, skip=False):
    """
    This is a stub that is used to keep track of what I need to still do.

    """
    substructures = load_protonation_substructs()
    out_containers = []
    inputs = [(cont, substructures) for cont in self.contnrs]

    tmp = mp.MultiThreading(inputs, self.params["num_processors"], parallel_addH)

    tmp = mp.flatten_list(tmp)
    contnr_indx_no_touch = Utils.contnrs_no_touchd(
        self, tmp
    )

    for miss_indx in contnr_indx_no_touch:
        Utils.log(
            "\tWARNING: OBABEL produced no valid protonation states for " +
            self.contnrs[miss_indx].orig_smi + " (" +
            self.contnrs[miss_indx].name + "), so using the original " +
            "smiles."
        )
        amol = self.contnrs[miss_indx].mol_orig_smi
        amol.contnr_idx = miss_indx

        amol.genealogy = [
            amol.orig_smi + " (source)",
            amol.orig_smi_deslt + " (desalted)",
            "(WARNING: OBABEL could not assign protonation states)"
        ]

        tmp.append(amol)

    ChemUtils.bst_for_each_contnr_no_opt(self, tmp)



def parallel_addH(container, substructures):
    """
    We take a container and a list of substructures and return all the appropriate protonation
    variants.

    :params container container: A container for a
    """
    return_value = []

    my_mol = container.mols[0]

    orig_mol = Chem.AddHs(Chem.RemoveHs(my_mol.rdkit_mol))
    mols = [orig_mol]
    protonation_sites = get_protonation_sites(orig_mol, substructures)

    for site in protonation_sites:
        mols = protonate_site(mols, site)

    smis = [Chem.MolToSmiles(mol) for mol in mols]
    rdkit_mols = [Chem.MolFromSmiles(smi) for smi in smis]

    # Convert from rdkit mols to MyMols and remove those with odd substructures
    addH_mols = [MyMol.MyMol(mol) for mol in rdkit_mols if mol is not None]
    addH_mols = [mol for mol in addH_mols if mol.crzy_substruc() == False]

    # I once saw it add a C+ here. So do a sanity check at
    # this point.
    for Hm in addH_mols:
        Hm.inherit_contnr_props(container)
        Hm.genealogy = my_mol.genealogy[:]
        Hm.name = my_mol.name

        if Hm.smiles() != my_mol.smiles():
            Hm.genealogy.append(Hm.smiles(True) + " (protonated)")

        return_value.append(Hm)

    return return_value
