
import sys

import gypsum.Multiprocess as mp
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils

from gypsum.MyMol import MyConformer

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)

try:
    import numpy
except:
    Utils.log("You need to install numpy and its dependencies.")
    sys.exit(0)

try:
    from scipy.cluster.vq import kmeans2
except:
    Utils.log("You need to install scipy and its dependencies.")
    sys.exit(0)


def GetRingConfs(mol, thoroughness, max_variants_per_compound):
    contnr_idx = mol.contnr_idx

    # All the ones in this contnr must have nonatomatic rings.
    # So just make a new mols list.

    # Get the ring atom indecies
    rings = mol.m_num_nonaro_rngs()

    # Convert that into the bond indecies.
    rings_by_bonds = []
    for ring_atom_indecies in rings:
        bond_indecies = []
        for ring_atm_idx in ring_atom_indecies:
            a = mol.GetAtomWithIdx(ring_atm_idx)
            bonds = a.GetBonds()
            for bond in bonds:
                atom_indecies = [
                    bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                ]
                atom_indecies.remove(ring_atm_idx)
                other_atm_idx = atom_indecies[0]
                if other_atm_idx in ring_atom_indecies:
                    bond_indecies.append(bond.GetIdx())
        bond_indecies = list(set(bond_indecies))
        bond_indecies.sort()

        rings_by_bonds.append(bond_indecies)

    # Generate a bunch of conformations, ordered from best
    # energy to worst. Note that this is cached.
    # Minimizing too.
    mol.add_conformers(
        thoroughness * max_variants_per_compound,
        0.1, True
    )

    if mol.GetNumConformers() > 0:
        # Sometimes there are no conformers if it's an impossible structure.
        # Like [H]c1nc(N2C(=O)[C@@]3(C([H])([H])[H])[C@@]4([H])O[C@@]([H])(C([H])([H])C4([H])[H])[C@]3(C([H])([H])[H])C2=O)sc1[H]
        # So don't save this one anyway.

        # Get the scores (lowest energy) of these minimized conformers
        mol.load_conformations_into_mol_3d()

        # Extract just the rings.
        ring_mols = [Chem.PathToSubmol(mol.rdkit_mol, bi)
                        for bi in rings_by_bonds]

        # Align get the rmsds relative to the first conformation,
        # for each ring separately.
        list_of_rmslists = [[]] * len(ring_mols)
        for k in range(len(ring_mols)):
            list_of_rmslists[k] = []
            AllChem.AlignMolConformers(
                ring_mols[k], RMSlist=list_of_rmslists[k]
            )

        # Get points for each conformer (rmsd_ring1, rmsd_ring2,
        # rmsd_ring3)
        pts = numpy.array(list_of_rmslists).T
        pts = numpy.vstack((numpy.array([[0.0] * pts.shape[1]]), pts))

        # cluster those points, get lowest-energy member of each
        if len(pts) < max_variants_per_compound:
            num_clusters = len(pts)
        else:
            num_clusters = max_variants_per_compound

        # try:
        groups = kmeans2(pts, num_clusters, minit='points')[1]

        #conformers = new_mymol
        # except:
        #     print pts
        #     print params["max_variants_per_compound"]
        #     print kmeans2(pts, num_clusters, minit='points')[1]

        # Note that you have some geometrically diverse conformations
        # here, but there could be other versions (enantiomers,
        # tautomers, etc.) that also contribute similar conformations.
        # In the end, you'll be selecting from all these together, so
        # similar ones could end up together.

        best_ones = {}
        conformers = mol.rdkit_mol.GetConformers()
        for k, grp in enumerate(groups):
            if not grp in best_ones.keys():
                best_ones[grp] = mol.conformers[k]
        best_confs = best_ones.values()

        # Convert rdkit mols to MyMol
        results = []
        for conf in best_confs:
            new_mol = mol.copy()
            c = MyConformer(new_mol, conf.conformer())
            new_mol.conformers = [c]
            energy = c.energy

            new_mol.genealogy = mol.genealogy[:]
            new_mol.genealogy.append(
                new_mol.smiles(True) +
                " (nonaromatic ring conformer: " + str(energy) +
                " kcal/mol)"
            )

            results.append(new_mol)  # i is mol index

        return results



def generate_alternate_3d_nonaromatic_ring_confs(contnrs, thoroughness, max_variants_per_compound, num_processors):
    """
    Docking programs like Vina rotate chemical moieties around their rotatable
    bonds, so it's not necessary to generate a larger rotomer library for each
    molecule. The one exception to this rule is non-aromatic rings, which can
    assume multiple conformations (boat vs. chair, etc.). This function
    generates a few low-energy ring structures for each molecule with a
    non-aromatic ring(s).
    """

    Utils.log(
        "Generating several conformers of molecules with non-aromatic " +
        "rings (boat vs. chair, etc.)..."
    )

    # params = []
    # for contnr in self.contnrs:
    #     for mol in contnr.mols:
    #         params.append((mol, self.params))

    params = []
    ones_with_nonaro_rngs = set([])
    for contnr_idx, contnr in enumerate(contnrs):
        if contnr.num_nonaro_rngs > 0:
            ones_with_nonaro_rngs.add(contnr_idx)
            for mol in contnr.mols:
                params.append((mol, thoroughness, max_variants_per_compound))

    if len(ones_with_nonaro_rngs) == 0:
        return  # There are no such ligands to process.

    #Utils.log("\tApplies to molecule derived from " + orig_smi)
    tmp = mp.MultiThreading(
        params, num_processors, GetRingConfs
    )

    results = [item for sublist in tmp for item in sublist]

    # # Remove mol list for the ones with nonaromatic rings
    # for contnr_idx in ones_with_nonaro_rngs:
    #     self.contnrs[contnr_idx].mols = []

    # Group by mol. You can't use existing functions because they would
    # require you to recalculate already calculated energies.
    grouped = {}
    for mol in results:
        # Save the energy as a prop while you're here.
        energy = mol.conformers[0].energy
        mol.mol_props["Energy"] = energy

        contnr_idx = mol.contnr_idx
        if not contnr_idx in grouped:
            grouped[contnr_idx] = []
        grouped[contnr_idx].append((energy, mol))

    # Now, for each contnr, keep only the best ones.
    for contnr_idx in grouped:
        lst = grouped[contnr_idx]

        if len(lst) != 0:
            contnrs[contnr_idx].mols = []
            lst.sort()
            lst = lst[:max_variants_per_compound]
            for energy, mol in lst:
                contnrs[contnr_idx].add_mol(mol)
        else:
            for i in range(len(contnrs[contnr_idx].mols)):
                contnrs[contnr_idx].mols[i].genealogy.append(
                    "(WARNING: Could not generate alternate conformations " +
                    "of nonaromatic ring)"
                )


