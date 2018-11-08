"""
This module is made to identify and enumerate the possible protonation sites of molecules.
"""

from rdkit import Chem

import gypsum.parallelizer as parallelizer
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol
import gypsum.MolContainer as MolCont

from gypsum.Steps.SMILES.protonation.protonation_functions import protonate

def add_hydrogens(contnrs, min_pH, max_pH, st_dev, max_variants,
                thoroughness, num_processors, multithread_mode, Parallelizer_obj):
    """
    This is a stub that is used to keep track of what I need to still do.

    """

    protonation_settings = {"min_ph": min_pH,
                            "max_ph": max_pH,
                            "st_dev": st_dev}

    inputs = [[cont, protonation_settings] for cont in contnrs]

    tmp = Parallelizer_obj.run(parallel_addH, inputs, num_processors, multithread_mode)

    tmp = parallelizer.flatten_list(tmp)

    contnr_indx_no_touch = Utils.contnrs_no_touchd(contnrs, tmp)

    for miss_indx in contnr_indx_no_touch:
        Utils.log(
            "\tWARNING: Gypsum produced no valid protonation states for " +
            contnrs[miss_indx].orig_smi + " (" +
            contnrs[miss_indx].name + "), so using the original " +
            "smiles."
        )
        amol = contnrs[miss_indx].mol_orig_smi
        amol.contnr_idx = miss_indx

        amol.genealogy = [
            amol.orig_smi + " (source)",
            amol.orig_smi_deslt + " (desalted)",
            "(WARNING: Gypsum could not assign protonation states)"
        ]

        tmp.append(amol)

    ChemUtils.bst_for_each_contnr_no_opt(contnrs, tmp, max_variants, thoroughness)

def parallel_addH(container, protonation_settings):
    """
    We take a container and a list of substructures and return all the
    appropriate protonation variants.

    :params container container: A container for a
    """
    return_value = []

    protonation_settings["smiles"] = container.orig_smi_canonical
    smis = protonate(protonation_settings)
    rdkit_mols = [Chem.MolFromSmiles(smi.strip()) for smi in smis]

    # Convert from rdkit mols to MyMols and remove those with odd substructures
    addH_mols = [MyMol.MyMol(mol) for mol in rdkit_mols if mol is not None]
    addH_mols = [mol for mol in addH_mols if mol.remove_bizarre_substruc() is False]

    # I once saw it add a C+ here. So do a sanity check at
    # this point.
    orig_mol = container.mol_orig_smi
    for Hm in addH_mols:
        Hm.inherit_contnr_props(container)
        Hm.genealogy = orig_mol.genealogy[:]
        Hm.name = orig_mol.name

        if Hm.smiles() != orig_mol.smiles():
            Hm.genealogy.append(Hm.smiles(True) + " (protonated)")

        return_value.append(Hm)

    return return_value
