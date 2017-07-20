import gypsum.Utils as Utils
import sys

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)

def save_to_sdf(self):
    """
    Saves the 3D models to the disk as an SDF file.
    """

    if self.params["separate_output_files"] == False:
        w = Chem.SDWriter(self.params["output_file"])
    else:
        w = Chem.SDWriter(self.params["output_file"] + ".params.sdf")

    # Save an empty molecule too with the parameters.
    m = Chem.Mol()
    m.SetProp("_Name", "EMPTY MOLECULE DESCRIBING GYPSUM PARAMETERS")
    for param in self.params:
        m.SetProp(param, str(self.params[param]))
    w.write(m)

    if self.params["separate_output_files"] == True:
        w.flush()
        w.close()

    Utils.log("Saving molecules associated with...")
    for i, contnr in enumerate(self.contnrs):
        contnr.add_container_properties()
        Utils.log("\t" + contnr.orig_smi)

        if self.params["separate_output_files"] == True:
            w = Chem.SDWriter(self.params["output_file"] + "." + str(i + 1) + ".sdf")

        for m in contnr.mols:
            m.load_conformations_into_mol_3d()
            w.write(m.rdkit_mol)

        if self.params["separate_output_files"] == True:
            w.flush()
            w.close()

    if self.params["separate_output_files"] == False:
        w.flush()
        w.close()

