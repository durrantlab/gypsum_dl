"""A module for generating alternate chiralities."""

import __future__

import copy
import itertools
import random

import gypsum.Parallelizer as Parallelizer
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

def enumerate_chiral_molecules(contnrs, max_variants_per_compound, thoroughness, num_procs, multithread_mode, parallelizer_obj):
    """Enumerates all possible enantiomers of a molecule. If the chirality of
       an atom is given, that chiral center is not varied. Only the chirality
       of unspecified chiral centers is varied.

    :param contnrs: A list of containers.
    :type contnrs: list
    :param max_variants_per_compound: [description] JDD: Figure out.
    :type max_variants_per_compound: int
    :param thoroughness: [description] JDD: Figure out.
    :type thoroughness: int
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param multithread_mode: The multithread mode.
    :type multithread_mode: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    # No point in continuing none requested.
    if max_variants_per_compound == 0:
        return

    Utils.log("Enumerating all possible enantiomers for all molecules...")

    # Group the molecules so you can feed them to parallelizer.
    params = []
    for contnr in contnrs:
        for mol in contnr.mols:
            params.append(tuple([mol, thoroughness, max_variants_per_compound]))
    params = tuple(params)

    # Run it through the parallelizer.
    tmp = parallelizer_obj.run(params, parallel_get_chiral, num_procs, multithread_mode)

    # Remove Nones (failed molecules)
    clean = Parallelizer.strip_none(tmp)

    # Flatten the data into a single list.
    flat = Parallelizer.flatten_list(clean)

    # Get the indexes of the ones that failed to generate.
    contnr_idxs_of_failed = Utils.find_missing_mol_idxs(contnrs, flat)

    # Go through the missing ones and throw a message.
    for miss_indx in contnr_idxs_of_failed:
        Utils.log(
            "\tCould not generate valid enantiomers for " +
            contnrs[miss_indx].orig_smi + " (" +
            contnrs[miss_indx].name + "), so using existing " +
            "(unprocessed) structures.")
        for mol in contnrs[miss_indx].mols:
            mol.genealogy.append("(WARNING: Unable to generate enantiomers)")
            clean.append(mol)

    # Keep only the top few compound variants in each container, to prevent a
    # combinatorial explosion.
    ChemUtils.bst_for_each_contnr_no_opt(contnrs, flat, max_variants_per_compound, thoroughness)

def parallel_get_chiral(mol, max_variants_per_compound, thoroughness):
    """A parallelizable function for enumerating chiralities.

    :param mol: The input molecule.
    :type mol: MyMol.MyMol
    :param max_variants_per_compound: [description] JDD: Figure out.
    :type max_variants_per_compound: int
    :param thoroughness: [description] JDD: Figure out.
    :type thoroughness: int
    :return: A list of MyMol.MyMol objects.
    :rtype: list
    """

    # Get all chiral centers that aren't assigned explicitly in the input
    # molecules.
    unasignd = [p[0] for p in mol.chiral_cntrs_w_unasignd() if p[1] == "?"]
    num = len(unasignd)

    # Get all possible chiral assignments. If the chirality is specified,
    # retain it.
    results = []
    if num == 0:
        # There are no unspecified chiral centers, so just keep existing.
        results.append(mol)
        return results
    elif num == 1:
        # There's only one chiral center.
        options = ["R", "S"]
    else:
        # There are multiple chiral centers.
        starting = [["R"], ["S"]]
        options = [["R"], ["S"]]
        for i in range(num - 1):
            options = list(itertools.product(options, starting))
            options = [list(itertools.chain(c[0], c[1])) for c in options]

    # Let the user know the number of chiral centers.
    Utils.log(
        "\t" + mol.smiles(True) + " (" + mol.name + ") has " +
        str(len(options)) + " enantiomers when chiral centers with " +
        "no specified chirality are systematically varied."
    )

    # Randomly select a few of the chiral combinations to examine. This is to
    # reduce the potential  combinatorial explosion.
    num_to_keep_initially = thoroughness * max_variants_per_compound
    options = Utils.random_sample(options, num_to_keep_initially, "")

    # Go through the chirality combinations and make a molecule with that
    # chirality.
    for option in options:
        # Copy the initial rdkit molecule≥
        a_rd_mol = copy.copy(mol.rdkit_mol)

        # Set its chirality.
        for idx, chiral in zip(unasignd, option):
            if chiral == "R":
                a_rd_mol.GetAtomWithIdx(idx).SetChiralTag(
                    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
                )
            elif chiral == "S":
                a_rd_mol.GetAtomWithIdx(idx).SetChiralTag(
                    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
                )

        # Make a new MyMol.MyMol object from that rdkit molecule.
        new_mol = MyMol.MyMol(a_rd_mol)

        # Add the new molecule to the list of results, if it does not have a
        # bizzare substructure.
        if not new_mol.remove_bizarre_substruc():
            new_mol.contnr_idx = mol.contnr_idx
            new_mol.name = mol.name
            new_mol.genealogy = mol.genealogy[:]
            new_mol.genealogy.append(
                new_mol.smiles(True) + " (chirality)"
            )
            results.append(new_mol)

    # Return the results.
    return results
