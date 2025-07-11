"""
This module identifies and enumerates the possible protonation sites of
molecules.
"""

from dimorphite_dl import protonate_smiles
from rdkit import Chem

import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl import MyMol, chem_utils, utils


def add_hydrogens(
    contnrs,
    min_pH,
    max_pH,
    st_dev,
    max_variants_per_compound,
    thoroughness,
    num_procs,
    job_manager,
    parallelizer_obj,
):
    """Adds hydrogen atoms to molecule containers, as appropriate for a given
       pH.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: A list.
    :param min_pH: The minimum pH to consider.
    :type min_pH: float
    :param max_pH: The maximum pH to consider.
    :type max_pH: float
    :param st_dev: The standard deviation. See Dimorphite-DL paper.
    :type st_dev: float
    :param max_variants_per_compound: To control the combinatorial explosion,
       only this number of variants (molecules) will be advanced to the next
       step.
    :type max_variants_per_compound: int
    :param thoroughness: How many molecules to generate per variant (molecule)
       retained, for evaluation. For example, perhaps you want to advance five
       molecules (max_variants_per_compound = 5). You could just generate five
       and advance them all. Or you could generate ten and advance the best
       five (so thoroughness = 2). Using thoroughness > 1 increases the
       computational expense, but it also increases the chances of finding good
       molecules.
    :type thoroughness: int
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param job_manager: The multithred mode to use.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    utils.log("Ionizing all molecules...")

    # Make a simple directory with the ionization parameters.
    protonation_settings = {
        "ph_min": min_pH,
        "ph_max": max_pH,
        "precision": st_dev,
        "max_variants": thoroughness * max_variants_per_compound,
    }

    # Format the inputs for use in the parallelizer.
    inputs = tuple(
        (cont, protonation_settings)
        for cont in contnrs
        if isinstance(cont.orig_smi_canonical, str)
    )

    # Run the parallelizer and collect the results.
    results = []
    if parallelizer_obj is None:
        results.extend(parallel_add_H(i[0], i[1]) for i in inputs)
    else:
        results = parallelizer_obj.run(inputs, parallel_add_H, num_procs, job_manager)
    results = Parallelizer.flatten_list(results)

    # Dimorphite-DL might not have generated ionization states for some
    # molecules. Identify those that are missing.
    contnr_idxs_of_failed = utils.fnd_contnrs_not_represntd(contnrs, results)

    # For those molecules, just use the original SMILES string, with hydrogen
    # atoms added using RDKit.
    for miss_indx in contnr_idxs_of_failed:
        utils.log(
            "\tWARNING: Gypsum-DL produced no valid ionization states for "
            + contnrs[miss_indx].orig_smi
            + " ("
            + contnrs[miss_indx].name
            + "), so using the original "
            + "smiles."
        )

        amol = contnrs[miss_indx].mol_orig_frm_inp_smi
        amol.contnr_idx = miss_indx

        # Save this failure to the genealogy record.
        amol.genealogy = [
            f"{amol.orig_smi} (source)",
            f"{amol.orig_smi_deslt} (desalted)",
            "(WARNING: Gypsum-DL could not assign ionization states)",
        ]

        # Save this one to the results too, even though not processed
        # properly.
        results.append(amol)

    # Keep only the top few compound variants in each container, to prevent a
    # combinatorial explosion.
    chem_utils.bst_for_each_contnr_no_opt(
        contnrs, results, max_variants_per_compound, thoroughness
    )


def parallel_add_H(contnr, protonation_settings):
    """Creates alternate ionization variants for a given molecule container.
       This is the function that gets fed into the parallelizer.

    :param contnr: The molecule container.
    :type contnr: MolContainer.MolContainer
    :param protonation_settings: Protonation settings to pass to Dimorphite-DL.
    :type protonation_settings: dict
    :return: [description]
    :rtype: [type]
    """

    # Make sure the canonical SMILES is actually a string.
    if not isinstance(contnr.orig_smi_canonical, str):
        utils.log(f"container.orig_smi_canonical: {contnr.orig_smi_canonical}")
        utils.log(
            f"type container.orig_smi_canonical: {str(type(contnr.orig_smi_canonical))}"
        )
        utils.exception(f"container.orig_smi_canonical: {contnr.orig_smi_canonical}")

    # Add the SMILES string to the protonation parameters.
    protonation_settings["smiles_input"] = contnr.orig_smi_canonical

    # Protonate the SMILESstring. This is Dimorphite-DL.
    smis = protonate_smiles(**protonation_settings)

    # Convert the protonated SMILES strings into a list of rdkit molecule
    # objects.
    rdkit_mols = [Chem.MolFromSmiles(smi.strip()) for smi in smis]

    # Convert from rdkit mols to MyMol.MyMol.
    addH_mols = [MyMol.MyMol(mol) for mol in rdkit_mols if mol is not None]

    # Remove MyMols with odd substructures.
    addH_mols = [mol for mol in addH_mols if mol.remove_bizarre_substruc() is False]

    # I once saw it add a "C+"" here. So do a secondary check at this point to
    # make sure it's valid. Recreate the list, moving new MyMol.MyMol objects
    # into the return_values list.

    return_values = []

    orig_mol = contnr.mol_orig_frm_inp_smi
    for Hm in addH_mols:
        Hm.inherit_contnr_props(contnr)
        Hm.genealogy = orig_mol.genealogy[:]
        Hm.name = orig_mol.name

        if Hm.smiles() != orig_mol.smiles():
            Hm.genealogy.append(f"{Hm.smiles(True)} (protonated)")

        return_values.append(Hm)

    return return_values
