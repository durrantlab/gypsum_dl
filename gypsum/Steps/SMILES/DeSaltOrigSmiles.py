from ... import ChemUtils
from ... import Utils
from ... import multiprocess_v2 as mp
from ... import MyMol

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)

def desalt_orig_smi(self):
    """
    If an input molecule has multiple unconnected fragments, this removes all
    but the largest fragment.
    """

    Utils.log(
        "Desalting all molecules (i.e., keeping only largest fragment)."
    )

    params = self.contnrs

    class DeSalter(mp.GeneralTask):
        def value_func(self, contnr, results_queue):
            # Split it into fragments
            frags = contnr.get_frags_of_orig_smi()

            if len(frags) == 1:
                # It's only got one fragment, so default assumption that
                # orig_smi = orig_smi_deslt is correct.
                return
            else:
                Utils.log(
                    "\tMultiple fragments found in " + contnr.orig_smi + 
                    " (" + contnr.name + ")"
                )
                
                # Find the biggest fragment
                num_heavy_atoms = []
                num_heavy_atoms_to_frag = {}

                for i, f in enumerate(frags):
                    num = f.GetNumHeavyAtoms()
                    num_heavy_atoms.append(num)
                    num_heavy_atoms_to_frag[num] = f

                max_num = max(num_heavy_atoms)
                biggest_frag = num_heavy_atoms_to_frag[max_num]

                # Return info about that biggest fragment.
                new_mol = MyMol.MyMol(biggest_frag)
                new_mol.contnr_idx = contnr.contnr_idx
                new_mol.name = contnr.name
                new_mol.makeMolFromSmiles() # Need to update the mol.

                self.results.append(new_mol)
    tmp = mp.MultiThreading(params, self.params["num_processors"], DeSalter)
    
    # Go through each contnr and update the orig_smi
    for desalt_mol in tmp.results:
        if self.contnrs[desalt_mol.contnr_idx].orig_smi != desalt_mol.orig_smi:
            self.contnrs[desalt_mol.contnr_idx].update_orig_smi(
                desalt_mol.orig_smi
            )



