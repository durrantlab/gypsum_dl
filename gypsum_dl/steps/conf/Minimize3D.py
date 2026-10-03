"""
This module performs a final 3D minimization to improve the small-molecule
geometry.
"""

import __future__

import copy
import operator

import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl import chem_utils, utils
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.MyMol import MyConformer


def minimize_3d(
    contnrs: list[MolContainer],
    max_variants_per_compound: int,
    thoroughness: int,
    num_procs: int,
    second_embed: bool,
    job_manager: str,
    parallelizer_obj: "Parallelizer.Parallelizer | None",
    include_nonaro_rings: bool = False,
    ring_conf_failed_contnr_idxs: frozenset[int] = frozenset(),
) -> None:
    """This function minimizes a 3D molecular conformation. In an attempt to
       not get trapped in a local minimum, it actually generates a number of
       conformers, minimizes them all, and then saves the best one.

    :param contnrs: A list of containers (MolContainer.MolContainer).
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
    :param job_manager: The multithred mode to use.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    :param include_nonaro_rings: Whether to also minimize molecules that have
        non-aromatic rings. Those molecules are minimized as a side effect of
        generating alternate ring conformations, so they only need to be
        handled here when that step did not run. Defaults to False.
    :type include_nonaro_rings: bool
    :param ring_conf_failed_contnr_idxs: Containers with non-aromatic rings
        that the ring-conformer step could not process, so it never minimized
        them. These are minimized here even when include_nonaro_rings is
        False. Defaults to an empty set.
    :type ring_conf_failed_contnr_idxs: frozenset[int]
    """

    # Let the user know you're on this step.
    utils.log("Minimizing all 3D molecular structures...")

    # Create the parameters (inputs) for the parallelizer.
    params = []
    ones_without_nonaro_rngs = set([])
    for contnr in contnrs:
        if (
            include_nonaro_rings
            or contnr.num_nonaro_rngs == 0
            or contnr.contnr_idx in ring_conf_failed_contnr_idxs
        ):
            # Unless the caller asks for them, ones with nonaromatic rings are
            # skipped here because generating their alternate ring
            # conformations already minimized them.
            for mol in contnr.mols:
                ones_without_nonaro_rngs.add(mol.contnr_idx)
                params.append(
                    (mol, max_variants_per_compound, thoroughness, second_embed)
                )
    params = tuple(params)

    # Run the inputs through the parallelizer.
    tmp = []
    if parallelizer_obj is None:
        # MultiThreading with one processor is the same dispatch serial mode
        # uses, and it draws a seed per job; calling run_one directly left
        # this path drawing from whatever generator state happened to be in
        # place, so --random_seed did not reach it.
        tmp = Parallelizer.MultiThreading(params, 1, parallel_minit)
    else:
        tmp = parallelizer_obj.run(params, parallel_minit, num_procs, job_manager)

    # Save energy into MyMol object, and get a list of just those objects.
    contnr_list_not_empty = set([])  # To keep track of which container lists
    # are not empty. These are the ones
    # you'll be repopulating with better
    # optimized structures.
    results = []  # Will contain MyMol.MyMol objects, with the saved energies
    # inside.
    for mol in Parallelizer.strip_none(tmp):
        mol.mol_props["Energy"] = utils.energy_for_output(mol.conformers[0].energy)
        results.append(mol)
        contnr_list_not_empty.add(mol.contnr_idx)

    # Go through each of the containers that are not empty and remove current
    # ones. Because you'll be replacing them with optimized versions.
    contnr_by_idx = utils.contnrs_by_idx(contnrs)
    for i in contnr_list_not_empty:
        contnr_by_idx[i].mols = []

    # Go through each of the minimized mols, and populate containers they
    # belong to.
    for mol in results:
        contnr_by_idx[mol.contnr_idx].add_mol(mol)

    # Alert the user to any errors, and drop the molecules behind them. Such a
    # molecule has nothing writable in it (load_conformers_into_rdkit_mol
    # returns early when rdkit_mol is None), but the steps that follow do call
    # accessors on it, in the main process, where an exception ends the run
    # after all the expensive work is done. Dropping it also lets
    # deal_with_failed_molecules report the input as failed.
    for contnr in contnrs:
        kept = []
        for mol in contnr.mols:
            if mol.rdkit_mol is None:
                mol.genealogy.append("(WARNING: Could not optimize 3D geometry)")
                mol.conformers = []
                continue
            kept.append(mol)
        contnr.mols = kept


def parallel_minit(mol, max_variants_per_compound, thoroughness, second_embed):
    """Minimizes the geometries of a MyMol.MyMol object. Meant to be run
    within parallelizer.

    :param mol: The molecule to minimize.
    :type mol: MyMol.MyMol
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
    :return: A molecule with the minimized conformers inside it.
    :rtype: MyMol.MyMol
    """

    # Serial and in-process runs pass the container's own molecule, while
    # multiprocessing passes a pickled copy. Work on a copy in every mode:
    # when every variant of a container fails here, minimize_3d keeps the
    # originals, and those must not have gained this attempt's conformers in
    # some job managers only.
    mol = copy.deepcopy(mol)

    # A cap of zero means "do not enumerate variants," not "emit no models,"
    # so one conformer is still needed here; a cap of zero would otherwise
    # request zero conformers and leave nothing to minimize.
    variant_cap = max(1, max_variants_per_compound)

    # Not minimizing. Just adding the conformers.
    mol.add_conformers(thoroughness * variant_cap, 0.1, False, second_embed)

    if len(mol.conformers) > 0:
        # Because it is possible to find a molecule that has no
        # acceptable conformers (i.e., is not possible geometrically).
        # Consider this:
        # O=C([C@@]1([C@@H]2O[C@@H]([C@@]1(C3=O)C)CC2)C)N3c4sccn4

        # Minimize every conformer rather than only the best few by
        # pre-minimization energy. That energy is dominated by embedding strain
        # and does not predict the minimized ranking (for flexible molecules it
        # is often anti-correlated), so pre-filtering on it discards good
        # conformers and leaves thoroughness with almost no effect when
        # max_variants_per_compound is small.
        minimized = mol.conformers
        for conf in minimized:
            conf.minimize()
        minimized.sort(key=operator.attrgetter("energy"))

        # Remove similar conformers
        # mol.eliminate_structurally_similar_conformers()

        # Get the best scoring (lowest energy) of these minimized conformers
        new_mol = copy.deepcopy(mol)
        c = MyConformer(new_mol, minimized[0].conformer(), second_embed)
        new_mol.conformers = [c]
        best_energy = c.energy

        # Save to the genealogy record.
        new_mol.genealogy = mol.genealogy[:]
        # Formatted rather than concatenated: smiles() reports failure as
        # None, and a TypeError here would be caught by the worker wrapper and
        # reported as a molecule that simply produced nothing.
        new_mol.genealogy.append(
            f"{new_mol.smiles(True)} (optimized conformer: "
            f"{utils.describe_energy(best_energy)})"
        )

        # Save best conformation. For some reason molecular properties
        # attached to mol are lost when returning from multiple
        # processors. So save the separately so they can be readded to
        # the molecule in a bit.
        # JDD: Still any issue?

        return new_mol
