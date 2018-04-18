"""Runs the smile preparation process."""
import __future__

from gypsum.Steps.SMILES.DeSaltOrigSmiles import desalt_orig_smi
from gypsum.Steps.SMILES.AddHydrogens import add_hydrogens
from gypsum.Steps.SMILES.MakeTautomers import make_tauts
from gypsum.Steps.SMILES.EnumerateChiralMols import enumerate_chiral_molecules
from gypsum.Steps.SMILES.EnumerateDoubleBonds import enumerate_double_bonds

def prepare_smiles(contnrs, params):
    """
    Runs the appropriate steps for processing the smile strings.
    """    
    
    print("Begin Desaltings")   
    max_variants_per_compound = params["max_variants_per_compound"]
    thoroughness = params["thoroughness"]
    num_processors = params["num_processors"]

    # Run the functions
    desalt_orig_smi(contnrs, num_processors)

    print("Done with Desalting")

    # Pausing for a moment to get that properly installed
    if not params["skip_adding_hydrogen"]:
        #JAKE@ REMOVED ADDHYDROGEN BECAUSE ITS BROKEN
        print("SKIPPING add_hydrogens(contnrs)")
    else:
        wrap_molecules(contnrs)

    if not params["skip_making_tautomers"]:
        make_tauts(contnrs, max_variants_per_compound, thoroughness, num_processors)
    
    if not params["skip_ennumerate_chiral_mol"]:
        enumerate_chiral_molecules(contnrs, max_variants_per_compound, thoroughness, num_processors)
    
    if not params["skip_ennumerate_double_bonds"]:
        enumerate_double_bonds(contnrs, max_variants_per_compound, thoroughness, num_processors)
    
def wrap_molecules(contnrs):
    """
    Problem: Each molecule container holds one smiles string
    (corresponding to the input structure). obabel produces multiple
    smiles strings at different pH values in the previous step. There is
    no way to store muliple smiles in a molecule container. But those
    containers are designed to store multiple RDKit molecule objects. To
    the previous step stores the differently protonated models as those
    objects, in the container's mol list.

    But, if the user skips the previous step, then the one smiles needs
    to be converted to a RDKit mol object for subsequent steps to work.
    Let's do that here.
    """
    for mol_cont in contnrs:
        if len(mol_cont.mols) == 0:
            smi = mol_cont.orig_smi_canonical
            mol_cont.add_smiles(smi)

