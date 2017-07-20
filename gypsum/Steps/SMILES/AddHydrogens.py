from ... import Multiprocess as mp
from ... import Utils
from ... import ChemUtils
from ... import MyMol
import random
import sys
import os
import tempfile

try:
    from rdkit.Chem import AllChem
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)


def add_hydrogens(self, skip=False):
    """
    Adds hydrogen atoms to the molecules, as appropriate for pH's ranging over
    the user-specified values. Note that though not a class function, it still
    accepts self as a parameter. This is the class that is calling it.
    """

    Utils.log("Adding hydrogen atoms...")

    tmp = None

    if skip:
        tmp = []
    else:
        # Save all smiles to a temporary file
        # PATRICK - We need to call the temp file function
        #flnm = str(random.randrange(0, 1000000)) + ".tmp"
        #while os.path.exists(flnm):
        #    flnm = str(random.randrange(0, 1000000)) + ".tmp"

        flnm = tempfile.mkstemp()[1]
        with open(flnm, 'w') as f:
            f.writelines(
                [m.orig_smi_deslt + " " + m.name + "____" + str(i) + "____" +
                m.orig_smi + "____" + m.orig_smi_deslt + "\n" for i, m in
                enumerate(self.contnrs)]
            )

        # Change pH in increments of 0.5 from min to max
        params = []
        pH = self.params["min_ph"]
        while pH <= self.params["max_ph"]:
            params.append((pH, flnm, self.params["openbabel_executable"]))
            pH = pH + self.params["delta_ph_increment"]

        tmp = mp.MultiThreading(params, self.params["num_processors"], parallel_addH)

        os.unlink(flnm)
        # *?* Check that this makes sense in the structure of ourput
        # Flattening return values
        tmp = mp.flatten_list(tmp)

    # Add back in remaining, using original smiles (no protonation). This is
    # better than nothing.
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


def parallel_addH(pH, flnm, obabel_loc):
    """
    A parallelizable helper function that adds hydrogens to molecules.

    :param float pH: The pH to consider when adding hydrogens
    :param str flnm: The file name that holds the file to consider.
    :param str obabel_loc: The location of the obabel installation.

    :results: Returns either the protonated RDKit molecule or a None object.
    """
    Utils.log("\tat pH " + str(pH))

    results = Utils.runit(
        obabel_loc + ' -p ' + str(pH) + ' -ismi ' +  flnm + ' -ocan'
    )

    return_value = []

    for s in results:
        s = s.strip()
        if s != "":
            prts = s.split()
            smi = prts[0]
            mol_info = " ".join(prts[1:])
            name, contnr_idx, orig_smi, orig_smi_deslt = mol_info.split("____")

            smi = fix_common_babel_ph_smiles_errors(smi)

            amol = MyMol.MyMol(smi)

            # I once saw it add a C+ here. So do a sanity check at
            # this point.
            if amol.rdkit_mol is not None:

                # Unfortuantely, obabel makes some systematic mistakes
                # when it comes to pH assignments. Try to correct them
                # here.
                # amol.fix_common_errors()

                if amol.crzy_substruc() == False:
                    # So no crazy substructure.

                    amol.contnr_idx = int(
                        contnr_idx
                    )

                    amol.genealogy.append(orig_smi + " (source)")

                    if orig_smi != orig_smi_deslt:
                        amol.genealogy.append(
                            orig_smi_deslt + " (desalted)"
                        )

                    amol.genealogy.append(
                        amol.smiles(True) + " (at pH " +
                        str(pH) + ")"
                    )

                    amol.name = name

                    return_value.append(amol)
                else:
                    # It has a crazy substructure.
                    Utils.log(
                        "\WARNING: " + smi + " (" + name  +
                        ") discarded."
                    )
    return return_value

def fix_common_babel_ph_smiles_errors(smiles):
    """
    Try to fix common structural erors.
    """
        
    mol = Chem.MolFromSmiles(smiles)

    while True:
        if mol is None:
            return None  # It's invalid somehow

        can_smi = Chem.MolToSmiles(
            mol, isomericSmiles=True, canonical=True
        )

        # inappropriate carbocations
        if "[C+]" in can_smi:
            mol = AllChem.ReplaceSubstructs(
                mol,
                Chem.MolFromSmarts('[$([C+](=*)(-*)-*)]'),
                Chem.MolFromSmiles('C')
            )[0]
            continue
        
        # Inappropriate modifications to carboxylic acids
        smrts = Chem.MolFromSmarts("C([O-])O")
        if mol.HasSubstructMatch(smrts):
            rxn = AllChem.ReactionFromSmarts(
                '[CH1:1](-[OH1:2])-[OX1-:3]>>[C:1](=[O:2])[O-:3]'
            )
            r = rxn.RunReactants([mol])
            if len(r) > 0:
                mol = r[0][0]
            continue

        # N+ bonded to only three atoms does not have a positive
        # charge.
        smrts = Chem.MolFromSmarts("[NX3+]")
        if mol.HasSubstructMatch(smrts):
            rxn = AllChem.ReactionFromSmarts('[NX3+:1]>>[N:1]')   
            r = rxn.RunReactants([mol])
            if len(r) > 0:
                mol = r[0][0]
            continue

        # If you get here, no changes were made, so return what you've got
        return can_smi


