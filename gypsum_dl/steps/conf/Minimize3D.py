"""
This module performs a final 3D minimization to improve the small-molecule
geometry.
"""

import copy

from gypsum_dl import utils
from gypsum_dl.molecule import MyConformer


def minimize_3d(
    contnrs,
    max_variants_per_compound,
    thoroughness,
    num_procs,
    second_embed,
    job_manager,
    parallelizer_obj,
):
    """This function minimizes a 3D molecular conformation. In an attempt to
       not get trapped in a local minimum, it actually generates a number of
       conformers, minimizes the best ones, and then saves the best of the
       best.

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
    :param job_manager: The multithred mode to use.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    # Let the user know you're on this step.
    utils.log("Minimizing all 3D molecular structures...")

    # Create the parameters (inputs) for the parallelizer.
    params = []
    ones_without_nonaro_rngs = set([])
    for contnr in contnrs:
        if contnr.num_nonaro_rngs == 0:
            # Because ones with nonaromatic rings have already been minimized,
            # so they can be skipped here.
            for mol in contnr.mols:
                ones_without_nonaro_rngs.add(mol.contnr_idx)
                params.append(
                    (mol, max_variants_per_compound, thoroughness, second_embed)
                )
    params = tuple(params)

    # Run the inputs through the parallelizer.
    tmp = []
    if parallelizer_obj is None:
        tmp.extend(parallel_minit(i[0], i[1], i[2], i[3]) for i in params)
    else:
        tmp = parallelizer_obj.run(params, parallel_minit, num_procs, job_manager)

    # Save energy into Molecule object, and get a list of just those objects.
    contnr_list_not_empty = set([])  # To keep track of which container lists
    # are not empty. These are the ones
    # you'll be repopulating with better
    # optimized structures.
    results = []  # Will contain Molecule objects, with the saved energies
    # inside.
    for mol in tmp:
        mol.mol_props["Energy"] = mol.conformers[0].energy
        results.append(mol)
        contnr_list_not_empty.add(mol.contnr_idx)

    # Go through each of the containers that are not empty and remove current
    # ones. Because you'll be replacing them with optimized versions.
    for i in contnr_list_not_empty:
        contnrs[i].mols = []

    # Go through each of the minimized mols, and populate containers they
    # belong to.
    for mol in results:
        contnrs[mol.contnr_idx].add_mol(mol)

    # Alert the user to any errors.
    for contnr in contnrs:
        for mol in contnr.mols:
            if mol.rdkit_mol == "":
                mol.genealogy.append("(WARNING: Could not optimize 3D geometry)")
                mol.conformers = []


def parallel_minit(mol, max_variants_per_compound, thoroughness, second_embed):
    """Minimizes the geometries of a Molecule object. Meant to be run
    within parallelizer.

    :param mol: The molecule to minimize.
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
    :return: A molecule with the minimized conformers inside it.
    :rtype: Molecule
    """

    # Not minimizing. Just adding the conformers.
    mol.add_conformers(thoroughness * max_variants_per_compound, 0.1, False)

    if len(mol.conformers) > 0:
        # Because it is possible to find a molecule that has no
        # acceptable conformers (i.e., is not possible geometrically).
        # Consider this:
        # O=C([C@@]1([C@@H]2O[C@@H]([C@@]1(C3=O)C)CC2)C)N3c4sccn4

        # Further minimize the unoptimized conformers that were among the best
        # scoring.
        max_vars_per_cmpd = max_variants_per_compound
        for i in range(len(mol.conformers[:max_vars_per_cmpd])):
            mol.conformers[i].minimize()

        # Remove similar conformers
        # mol.eliminate_structurally_similar_conformers()

        # Get the best scoring (lowest energy) of these minimized conformers
        new_mol = copy.deepcopy(mol)
        c = MyConformer(new_mol, mol.conformers[0].conformer(), second_embed)
        new_mol.conformers = [c]
        best_energy = c.energy

        # Save to the genealogy record.
        new_mol.genealogy = mol.genealogy[:]
        new_mol.genealogy.append(
            new_mol.smiles(True)
            + " (optimized conformer: "
            + str(best_energy)
            + " kcal/mol)"
        )

        # Save best conformation. For some reason molecular properties
        # attached to mol are lost when returning from multiple
        # processors. So save the separately so they can be readded to
        # the molecule in a bit.
        # JDD: Still any issue?

        return new_mol
