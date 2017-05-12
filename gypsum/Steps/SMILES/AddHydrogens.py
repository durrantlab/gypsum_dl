import gypsum.Multiprocess as mp
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol
import random
import sys
import os
import tempfile

try:
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

            amol = MyMol.MyMol(smi)

            # I once saw it add a C+ here. So do a sanity check at
            # this point.
            if amol.rdkit_mol is not None:

                # Unfortuantely, obabel makes some systematic mistakes
                # when it comes to pH assignments. Try to correct them
                # here.

                amol.fix_common_errors()

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
