from ... import multiprocess_v2 as mp
from ... import Utils
from ... import ChemUtils
from ... import MyMol
import random
import os

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)


def add_hydrogens(self):
    """
    Adds hydrogen atoms to the molecules, as appropriate for pH's ranging over
    the user-specified values. Note that though not a class function, it still
    accepts self as a parameter. This is the class that is calling it.
    """
    
    Utils.log("Adding hydrogen atoms...")

    # Save all smiles to a temporary file
    flnm = str(random.randrange(0,1000000)) + ".tmp"
    while os.path.exists(flnm):
        flnm = str(random.randrange(0,1000000)) + ".tmp"
    
    f = open(flnm, 'w')
    f.writelines(
        [m.orig_smi_deslt + " " + m.name + "____" + str(i) + "____" +
        m.orig_smi + "____" + m.orig_smi_deslt + "\n" for i, m in
        enumerate(self.contnrs)]
    )
    f.close()

    # Change pH in increments of 0.5 from min to max
    params = []
    pH = self.params["min_ph"]
    while pH <= self.params["max_ph"]:
        params.append((pH, f.name, self.params["openbabel_executable"]))
        pH = pH + self.params["delta_ph_increment"]

    class addH(mp.GeneralTask):
        def value_func(self, items, results_queue):
            pH, flnm, obabel_loc = items
            Utils.log("\tat pH " + str(pH))
            results = Utils.runit(
                obabel_loc + ' -p ' + str(pH) + ' -ismi ' +  f.name + ' -ocan'
            )

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

                            self.results.append(amol)
                        else:
                            Utils.log(
                                "\WARNING: " + smi + " (" + name  + 
                                ") discarded."
                            )
    tmp = mp.MultiThreading(params, self.params["num_processors"], addH)
        
    # Add back in remaining, using original smiles (no protonation). This is
    # better than nothing.

    contnr_indx_no_touch = Utils.contnrs_no_touchd(
        self, tmp.results
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

        tmp.results.append(amol)

    os.unlink(f.name)

    ChemUtils.bst_for_each_contnr_no_opt(self, tmp.results)
