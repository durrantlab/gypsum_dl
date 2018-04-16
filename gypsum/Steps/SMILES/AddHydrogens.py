"""
This module is made to identify and enumerate the possible protonation sites of molecules.
"""

import copy

from rdkit import Chem
import gypsum.Multiprocess as mp
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol

def load_protonation_substructs(min_pH=6.4, max_pH=8.4, pKa_std_range=1):
    """
    A pre-calculated list of R-groups with protonation sites, with their likely pKa bins.
    """
    subs = []
    with open("site_structures.smarts", 'r') as substruct:
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

                #NEED TO DIVIDE THIS BY 3s
                pKa_ranges = [splits[i:i+3] for i in range(3, len(splits)-1, 3)]
                #print(splits, '\n', pKa_ranges)

                prot = []
                for pKa_range in pKa_ranges:
                    site = pKa_range[0]
                    mean = pKa_range[1]
                    std = pKa_range[2]

                    protonation_state = define_protonation_state(mean, std, min_pH, \
                        max_pH, pKa_std_range)

                    prot.append([site, protonation_state])

                #!!! This takes any explicit hydrogen as a site for protonation.
                # We may find that we need to revisit this, as Nitrogen causes issues
                # with RDKit's kekulization when it has 4 bonds in aromatic rings.
                sub["prot"] = prot
                subs.append(sub)
    return subs

def define_protonation_state(mean, std, min_pH, max_pH, pKa_std_range):
    """
    Updates the substructure definitions to include the protonation state based on the user-given
    pH range. The size of the pKa range is also based on the number of standard deviations to be
    considered by the user param.
    """
    min_pKa = mean - (std * pKa_std_range)
    max_pKa = mean + (std * pKa_std_range)

    # This needs to be reassigned, and 'ERROR' should never make it past the next set of checks.
    protonation_state = 'ERROR'

    if min_pKa <= max_pH and min_pH <= max_pKa:
        protonation_state = 'BOTH'
    elif mean > max_pH:
        protonation_state = 'PROTONATED'
    elif mean < min_pH:
        protonation_state = 'DEPROTONATED'

    # We are error handling here
    if protonation_state == 'ERROR':
        print("HORRIBLE NONSENSE HAS OCCURED. MEAN: ",mean,", MIN_PH: ",min_pH,", MAX_PH: ",max_pH)
    
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
        keep_flag = True
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            protected = atom.GetProp("_protected")
            if protected == "1":
                keep_flag = False
        if keep_flag:
            unprotected_matches.append(match)
    return unprotected_matches

def get_protonation_sites(mol, subs):
    """
    For a single molecule, find all possible matches in the protonation R-group list,
    subs. Items that are higher on the list will be matched first, to the exclusion of
    later items.
    Returns a list of protonation sites and their pKa bin. ('Acid', 'Neutral', or 'Base')
    """
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
                    proton = site[0]
                    category = site[1]
                    new_site = (match[proton], category)
                    protonation_sites.append(new_site)
                protect_molecule(mol, match)

    return protonation_sites

def protonate_site(smis, site):
    idx, charge = site

    protonate = []
    deprotonate = []

    if charge == "DEPROTONATED" or charge == "BOTH":
        deprotonate = copy.deepcopy(smis)
    if charge == "PROTONATED" or charge == "BOTH":
        protonate = copy.deepcopy(smis)

    output_smis = []

    for smi in deprotonate:
        mol = Chem.MolFromSmiles(smi)
        atom = mol.GetAtomWithIdx(idx)
        element = atom.GetAtomicNum()
        if element == 7:
            atom.SetFormalCharge(0)
        else:
            atom.SetFormalCharge(-1)
        out_smile = Chem.MolToSmiles(mol)
        output_smis.append(out_smile)

    for smi in protonate:
        mol = Chem.MolFromSmiles(smi)
        atom = mol.GetAtomWithIdx(idx)
        element = atom.GetAtomicNum()
        if element == 7:
            atom.SetFormalCharge(+1)
        else:
            atom.SetFormalCharge(0)
        out_smile = Chem.MolToSmiles(mol)
        output_smis.append(out_smile)

    return output_smis

def add_hydrogens(self, skip=False):
    """
    This is a stub that is used to keep track of what I need to still do.

    """
    min_pH = self.params["min_pH"]
    max_pH = self.params["max_pH"]
    pKa_std_range = self.params["pKa_std_range"]

    substructures = load_protonation_substructs(min_pH, max_pH, pKa_std_range)
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
    We take a container and a list of substructures and return all the
    appropriate protonation variants.

    :params container container: A container for a
    """
    return_value = []

    orig_mol = container.Get_Mol()
    smis = [container.Get_Cannonical_Smile()]
    protonation_sites = get_protonation_sites(orig_mol, substructures)

    for site in protonation_sites:
        smis = protonate_site(smis, site)

    rdkit_mols = [Chem.MolFromSmiles(smi) for smi in smis]

    # Convert from rdkit mols to MyMols and remove those with odd substructures
    addH_mols = [MyMol.MyMol(mol) for mol in rdkit_mols if mol is not None]
    addH_mols = [mol for mol in addH_mols if mol.crzy_substruc() == False]

    # I once saw it add a C+ here. So do a sanity check at
    # this point.
    for Hm in addH_mols:
        Hm.inherit_contnr_props(container)
        Hm.genealogy = orig_mol.genealogy[:]
        Hm.name = orig_mol.name

        if Hm.smiles() != orig_mol.smiles():
            Hm.genealogy.append(Hm.smiles(True) + " (protonated)")

        return_value.append(Hm)

    return return_value

