"""
This module generates alternate non-aromatic ring conformations on the fly,
since most modern docking programs (e.g., Vina) can't consider alternate ring
conformations.
"""

import copy
import warnings

import numpy
from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.cluster.vq import kmeans2

import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl.molecule import MyConformer


def generate_alternate_3d_nonaromatic_ring_confs(
    contnrs,
    max_variants_per_compound,
    thoroughness,
    num_procs,
    second_embed,
    job_manager,
    parallelizer_obj,
):
    """Docking programs like Vina rotate chemical moieties around their
       rotatable bonds, so it's not necessary to generate a larger rotomer
       library for each molecule. The one exception to this rule is
       non-aromatic rings, which can assume multiple conformations (boat vs.
       chair, etc.). This function generates a few low-energy ring structures
       for each molecule with a non-aromatic ring(s).

    :param contnrs: A list of containers (container.MoleculeContainer).
    :type contnrs: list
    :param max_variants_per_compound: To control the combinatorial explosion,
       only this number of variants (molecules) will be advanced to the next
       step.
    :type max_variants_per_compound: int
    :param thoroughness: How many molecules to generate per variant (molecule)
       retained, for evaluation. For example, perhaps you want to advance five
       molecules (max_variants_per_compound = 5). You could just generate five
       and advance them all. Or you could generate ten and advance the best
       five (so thoroughness = 2). Using thoroughness > 1 increases the
       computational expense, but it also increases the chances of finding good
       molecules.
    :type thoroughness: int
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param second_embed: Whether to try to generate 3D coordinates using an
        older algorithm if the better (default) algorithm fails. This can add
        run time, but sometimes converts certain molecules that would
        otherwise fail.
    :type second_embed: bool
    :param job_manager: The multiprocess mode.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    :return: Returns None if no ring conformers are generated
    :rtype: None
    """

    # Let the user know you've started this step.
    logger.info(
        "Generating several conformers of molecules with non-aromatic "
        + "rings (boat vs. chair, etc.)..."
    )

    # Create parameters (inputs) to feed to the parallelizer.
    params = []
    ones_with_nonaro_rngs = set([])  # This is just to keep track of which
    # ones have non-aromatic rings.
    for contnr_idx, contnr in enumerate(contnrs):
        if contnr.num_nonaro_rngs > 0:
            ones_with_nonaro_rngs.add(contnr_idx)
            params.extend(
                (mol, max_variants_per_compound, thoroughness, second_embed)
                for mol in contnr.mols
            )
    params = tuple(params)

    # If there are no compounds with non-aromatic rings, no need to continue.
    if not ones_with_nonaro_rngs:
        return  # There are no such ligands to process.

    # Run it through the parallelizer
    tmp = []
    if parallelizer_obj is None:
        tmp.extend(parallel_get_ring_confs(i[0], i[1], i[2], i[3]) for i in params)
    else:
        tmp = parallelizer_obj.run(
            params, parallel_get_ring_confs, num_procs, job_manager
        )
    # Flatten the results.
    results = Parallelizer.flatten_list(tmp)

    # Group by mol. You can't use existing functions because they would
    # require you to recalculate already calculated energies.
    grouped = {}  # Index will be container index. Value is list of
    # (energy, mol) pairs.
    for mol in results:
        # Save the energy as a prop while you're here.
        energy = mol.conformers[0].energy
        mol.mol_props["Energy"] = energy

        # Add the mol with it's energy to the appropriate entry in grouped.
        # Make that entry if needed.
        contnr_idx = mol.contnr_idx
        if contnr_idx not in grouped:
            grouped[contnr_idx] = []
        grouped[contnr_idx].append((energy, mol))

    # Now, for each container, keep only the best ones.
    for contnr_idx, lst_enrgy_mol_pairs in grouped.items():
        if len(lst_enrgy_mol_pairs) != 0:
            contnrs[contnr_idx].mols = []  # Note that only affects ones that
            # had non-aromatic rings.
            lst_enrgy_mol_pairs.sort()  # Sorting by energy (first item in
            # pair).

            # Keep only the top ones.
            lst_enrgy_mol_pairs = lst_enrgy_mol_pairs[:max_variants_per_compound]

            # Add the top ones to the container mol list.
            for energy, mol in lst_enrgy_mol_pairs:
                contnrs[contnr_idx].add_mol(mol)
        else:
            # There are no entries in the list. It apparently wasn't able to
            # generate any alternate conformers. Let the user know.
            for i in range(len(contnrs[contnr_idx].mols)):
                contnrs[contnr_idx].mols[i].genealogy.append(
                    "(WARNING: Could not generate alternate conformations "
                    + "of nonaromatic ring)"
                )


def parallel_get_ring_confs(mol, max_variants_per_compound, thoroughness, second_embed):
    """Gets alternate ring conformations. Meant to run with the parallelizer class.

    :param mol: The molecule to process (with non-aromatic ring(s)).
    :type mol: Molecule
    :param max_variants_per_compound: To control the combinatorial explosion,
       only this number of variants (molecules) will be advanced to the next
       step.
    :type max_variants_per_compound: int
    :param thoroughness: How many molecules to generate per variant (molecule)
       retained, for evaluation. For example, perhaps you want to advance five
       molecules (max_variants_per_compound = 5). You could just generate five
       and advance them all. Or you could generate ten and advance the best
       five (so thoroughness = 2). Using thoroughness > 1 increases the
       computational expense, but it also increases the chances of finding good
       molecules.
    :type thoroughness: int
    :param second_embed: Whether to try to generate 3D coordinates using an
        older algorithm if the better (default) algorithm fails. This can add
        run time, but sometimes converts certain molecules that would
        otherwise fail.
    :type second_embed: bool
    :return: A list of Molecule objects, with alternate ring conformations.
    :rtype: list
    """

    # Make it easier to access the container index.
    contnr_idx = mol.contnr_idx

    # All the molecules in this container must have nonatomatic rings (because
    # they are all variants of the same source molecule). So just make a new
    # mols list.

    # Get the ring atom indecies
    rings = mol.get_idxs_of_nonaro_rng_atms()

    # Convert that into the bond indecies.

    # A list of lists, where each inner list has the indexes of the bonds that
    # comprise a ring.
    rings_by_bond_indexes = []
    for ring_atom_indecies in rings:
        bond_indexes = []
        for ring_atm_idx in ring_atom_indecies:
            a = mol.rdkit_mol.GetAtomWithIdx(ring_atm_idx)
            bonds = a.GetBonds()
            for bond in bonds:
                atom_indecies = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
                atom_indecies.remove(ring_atm_idx)
                other_atm_idx = atom_indecies[0]
                if other_atm_idx in ring_atom_indecies:
                    bond_indexes.append(bond.GetIdx())
        bond_indexes = sorted(set(bond_indexes))
        rings_by_bond_indexes.append(bond_indexes)

    # Generate a bunch of conformations, ordered from best energy to worst.
    # Note that this is cached. Minimizing too.
    mol.add_conformers(thoroughness * max_variants_per_compound, 0.1, True)

    if len(mol.conformers) > 0:
        # Sometimes there are no conformers if it's an impossible structure.
        # Like
        # [H]c1nc(N2C(=O)[C@@]3(C([H])([H])[H])[C@@]4([H])O[C@@]([H])(C([H])([H])C4([H])[H])[C@]3(C([H])([H])[H])C2=O)sc1[H]
        # So don't save this one anyway.

        # Get the scores (lowest energy) of these minimized conformers.
        mol.load_conformers_into_rdkit_mol()

        # Extract just the rings.
        ring_mols = [
            Chem.PathToSubmol(mol.rdkit_mol, bi) for bi in rings_by_bond_indexes
        ]

        # Align get the rmsds relative to the first conformation, for each
        # ring separately.
        list_of_rmslists = [[]] * len(ring_mols)
        for k in range(len(ring_mols)):
            list_of_rmslists[k] = []
            AllChem.AlignMolConformers(ring_mols[k], RMSlist=list_of_rmslists[k])

        # Get points for each conformer (rmsd_ring1, rmsd_ring2, rmsd_ring3)
        pts = numpy.array(list_of_rmslists).T
        pts = numpy.vstack((numpy.array([[0.0] * pts.shape[1]]), pts))

        # Cluster those points, get lowest-energy member of each.
        if len(pts) < max_variants_per_compound:
            num_clusters = len(pts)
        else:
            num_clusters = max_variants_per_compound

        # When kmeans2 runs on insufficient clusters, it can sometimes throw an
        # error about empty clusters. This is not necessary to throw for the
        # user and so we have supressed it here.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            groups = kmeans2(pts, num_clusters, minit="points")[1]

        # Note that you have some geometrically diverse conformations here, but
        # there could be other versions (enantiomers, tautomers, etc.) that also
        # contribute similar conformations. In the end, you'll be selecting from
        # all these together, so similar ones could end up together.

        # Key is group id from kmeans (int). Values are the MyConformers
        # objects.
        best_conf_per_group = {}

        conformers = mol.rdkit_mol.GetConformers()
        for k, grp in enumerate(groups):
            if grp not in list(best_conf_per_group.keys()):
                best_conf_per_group[grp] = mol.conformers[k]
        # best_confs has the MyConformers objects.
        best_confs = best_conf_per_group.values()

        # Convert rdkit mols to Molecule and save those Molecule objects
        # for returning.
        results = []
        for conf in best_confs:
            new_mol = copy.deepcopy(mol)
            c = MyConformer(new_mol, conf.conformer(), second_embed)
            new_mol.conformers = [c]
            energy = c.energy

            new_mol.genealogy = mol.genealogy[:]
            new_mol.genealogy.append(
                new_mol.smiles(True)
                + " (nonaromatic ring conformer: "
                + str(energy)
                + " kcal/mol)"
            )

            results.append(new_mol)  # i is mol index

        return results

    # If you get here, something went wrong.
    return None
