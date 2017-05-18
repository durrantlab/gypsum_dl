"""Runs the smile preparation process."""
from gypsum.Steps.SMILES.DeSaltOrigSmiles import desalt_orig_smi
from gypsum.Steps.SMILES.AddHydrogens import add_hydrogens
from gypsum.Steps.SMILES.MakeTautomers import make_tauts
from gypsum.Steps.SMILES.EnumerateChiralMols import enumerate_chiral_molecules
from gypsum.Steps.SMILES.EnumerateDoubleBonds import enumerate_double_bonds

def prepare_smiles(self):
    """
    Runs the appropriate steps for processing the smile strings.
    """
    desalt_orig_smi(self)

    if self.params["break"] == "desalt":
        save_partial_output(self)
        raise Exception
        ##############################

    if not self.params["skip_adding_hydrogen"]:
        add_hydrogens(self)
    else:
        wrap_molecules(self)
    self.print_current_smiles()

    if self.params["break"] == "add_hydrogens":
        save_partial_output(self)
        raise Exception
        ##############################

    if not self.params["skip_making_tautomers"]:
        make_tauts(self)
    self.print_current_smiles()

    if self.params["break"] == "make_tauts":
        save_partial_output(self)
        raise Exception
        ##############################

    if not self.params["skip_ennumerate_chiral_mol"]:
        enumerate_chiral_molecules(self)
    self.print_current_smiles()

    if self.params["break"] == "make_chiral":
        save_partial_output(self)
        raise Exception
        ##############################

    if not self.params["skip_ennumerate_double_bonds"]:
        enumerate_double_bonds(self)
    self.print_current_smiles()

    if self.params["break"] == "make_bonds":
        save_partial_output(self)
        raise Exception
        ##############################

def wrap_molecules(self):
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
    for i, mol_cont in enumerate(self.contnrs):
        if len(mol_cont.mols) == 0:
            smi = mol_cont.orig_smi_canonical
            mol_cont.add_smiles(smi)

def save_partial_output(self):
    output_file = self.params["output_file"]
    with open(output_file, 'w') as output:
        for i, mol_cont in enumerate(self.contnrs):
            output.write(str(i) + "\tOrig\t" + str(mol_cont.orig_smi_canonical)+ "\n")
            output.write(str(i) + "\tDesalt\t" + str(mol_cont.orig_smi_deslt)+ "\n")
            for smile in mol_cont.all_smiles():
                output.write(str(i) + "\tNonCan\t" + str(smile)+ "\n")
