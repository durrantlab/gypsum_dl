from ... import mp_queue as mp
from ... import Utils
from ... import ChemUtils
from ... MyMol import MyConformer
import copy



def minit(mol, params):
    
    # Utils.log("\tMinimizing one of the structures generated for " +
    # orig_smi)
    # Not minimizing
    mol.add_conformers(
        params["thoroughness"] * params["max_variants_per_compound"],
        0.1, False
    )

    if mol.GetNumConformers() > 0:
        # Because it is possible to find a molecule that has no
        # acceptable conformers (i.e., is not possible geometrically).
        # Consider this:
        # O=C([C@@]1([C@@H]2O[C@@H]([C@@]1(C3=O)C)CC2)C)N3c4sccn4
    
        # Further minimize the unoptimized conformers that were among
        # the best scoring.
        max_vars_per_cmpd = params["max_variants_per_compound"]
        for i in range(len(mol.conformers[:max_vars_per_cmpd])):
            mol.conformers[i].minimize()

        # Remove similar conformers
        #mol.eliminate_structurally_similar_conformers()

        # Get the best scoring (lowest energy) of these minimized
        # conformers
        new_mol = mol.copy()
        c = MyConformer(new_mol, mol.conformers[0].conformer())
        new_mol.conformers = [c]
        best_energy = c.energy

        new_mol.genealogy = mol.genealogy[:]
        new_mol.genealogy.append(
            new_mol.smiles(True) + " (optimized conformer: " +
            str(best_energy) + " kcal/mol)"
        )
        
        # Save best conformation. For some reason molecular properties
        # attached to mol are lost when returning from multiple
        # processors. So save the separately so they can be readded to
        # the molecule in a bit. props =
        return new_mol

def minimize_3d(self):
    """
    This function minimizes a 3D molecular conformation. In an attempt to not
    get trapped in a local minimum, it actually generates a number of
    conformers, minimizes the best ones, and then saves the best of the best.
    """

    Utils.log("Minimizing all 3D molecular structures...")

    params = []
    ones_without_nonaro_rngs = set([])
    for contnr in self.contnrs:
        if contnr.num_nonaro_rngs == 0:
            # Because ones with nonaromatic rings have already been minimized.
            for mol in contnr.mols:
                ones_without_nonaro_rngs.add(mol.contnr_idx)
                params.append((mol, self.params))
    
    tmp = mp.MultiThreading(params, self.params["num_processors"], minit)
    
    results = []
    # Save energy into MyMol object, and get a list of just those objects.
    contnr_list_not_empty = set([])
    for mol in tmp:
        mol.mol_props["Energy"] = mol.conformers[0].energy
        results.append(mol)
        contnr_list_not_empty.add(mol.contnr_idx)
    
    # Go through each of the contnrs, remove current ones
    for i in contnr_list_not_empty:  # ones_without_nonaro_rngs:
        self.contnrs[i].mols = []
    
    # Go through each of the minimized mols, and populate containers that way.
    for mol in results:
        self.contnrs[mol.contnr_idx].add_mol(mol)
    
    # Set the makeMol3D() flags that aren't set.
    for contnr_idx, contnr in enumerate(self.contnrs):
        for mol_index, mol in enumerate(contnr.mols):
            if mol.rdkit_mol == "":
                mol.genealogy.append(
                    "(WARNING: Could not optimize 3D geometry)"
                )
                mol.conformers = []
