"""Runs the smile preparation process. Generates alternate ionization,
tautomeric, chiral forms, etc."""

import __future__

from gypsum import Utils
from gypsum.Steps.SMILES.DeSaltOrigSmiles import desalt_orig_smi
from gypsum.Steps.SMILES.AddHydrogens import add_hydrogens
from gypsum.Steps.SMILES.MakeTautomers import make_tauts
from gypsum.Steps.SMILES.EnumerateChiralMols import enumerate_chiral_molecules
from gypsum.Steps.SMILES.EnumerateDoubleBonds import enumerate_double_bonds

def prepare_smiles(contnrs, params):
    """Runs the appropriate steps for processing the SMILES strings.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param params: The user parameters.
    :type params: dict
    """

    # Unpack some of the parameter values.
    min_ph = params["min_ph"]
    max_ph = params["max_ph"]
    std_dev = params["ph_std_dev"]
    max_variants_per_compound = params["max_variants_per_compound"]
    thoroughness = params["thoroughness"]
    num_procs = params["num_processors"]
    multithread_mode = params["multithread_mode"]
    parallelizer_obj = params["Parallelizer"]

    debug = False

    # Desalt the molecules.
    print("Begin Desaltings")
    desalt_orig_smi(contnrs, num_procs, multithread_mode, parallelizer_obj)
    print("Done with Desalting")

    if debug: Utils.print_current_smiles(contnrs)

    # Add hydrogens for user-specified pH, if requested.
    if not params["skip_adding_hydrogen"]:
        print("Protonating Molecules")
        add_hydrogens(contnrs, min_ph, max_ph, std_dev, max_variants_per_compound,
                      thoroughness, num_procs, multithread_mode,
                      parallelizer_obj)
        print("Done with Protonating")
    else:
        print("Skipping Protonation")
        wrap_molecules(contnrs)

    if debug: Utils.print_current_smiles(contnrs)

    # Make alternate tautomeric forms, if requested.
    if not params["skip_making_tautomers"]:
        print("Tautomerizing Molecules")
        make_tauts(contnrs, max_variants_per_compound, thoroughness,
                   num_procs, multithread_mode, parallelizer_obj)
        print("Done with Tautomerization")
    else:
        print("Skipping Tautomerization")

    if debug: Utils.print_current_smiles(contnrs)

    # Make alternate chiral forms, if requested.
    if not params["skip_ennumerate_chiral_mol"]:
        print("Enumerating Chirality")
        enumerate_chiral_molecules(contnrs, max_variants_per_compound,
                                   thoroughness, num_procs,
                                   multithread_mode, parallelizer_obj)
        print("Done with Chirality Enumeration")
    else:
        print("Skipping Chirality Enumeration")

    if debug: Utils.print_current_smiles(contnrs)

    # Make alternate double-bond isomers, if requested.
    if not params["skip_ennumerate_double_bonds"]:
        print("Enumerating Double Bonds")
        enumerate_double_bonds(contnrs, max_variants_per_compound,
                               thoroughness, num_procs,
                               multithread_mode, parallelizer_obj)
        print("Done with Double Bond Enumeration")
    else:
        print("Skipping Double Bond Enumeration")

    if debug: Utils.print_current_smiles(contnrs)

def wrap_molecules(contnrs):
    """Each molecule container holds only one SMILES string
    (corresponding to the input structure). Dimorphite-DL can potentially
    produce multiple SMILES strings with different ionization states. There is
    no way to store muliple smiles in a molecule container.

    But those containers are designed to store multiple RDKit molecule
    objects. Gypsum stores Dimorphite-DL output molecules as RdKit molecule
    objects, in the container's mol list.

    But Gypsum gives the user the option of skipping the protonation step. In
    this case, the one SMILES needs to be converted to a RDKit mol object for
    subsequent steps to work. Let's do that here.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    """

    for mol_cont in contnrs:
        if len(mol_cont.mols) == 0:
            smi = mol_cont.orig_smi_canonical
            mol_cont.add_smiles(smi)
