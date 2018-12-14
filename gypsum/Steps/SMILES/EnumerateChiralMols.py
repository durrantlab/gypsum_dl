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


def GetChiral(mol, max_variants_per_compound, thoroughness):

    # Get all chiral centers that aren't assigned.
    unasignd = [p[0] for p in mol.chiral_cntrs_w_unasignd()
                if p[1] == "?"]
    num = len(unasignd)
    results = []

    # Get all possible chiral assignments. If chiral is specified,
    # retain it.
    if num == 0:
        # There are no unspecified chiral centers, so just keep
        # existing.
        results.append(mol)
        return results
    elif num == 1:
        options = ["R", "S"]
    else:
        starting = [["R"], ["S"]]
        options = [["R"], ["S"]]
        for i in range(num - 1):
            options = list(itertools.product(options, starting))
            options = [list(itertools.chain(c[0], c[1]))
                       for c in options]

    Utils.log(
        "\t" + mol.smiles(True) + " (" + mol.name + ") has " +
        str(len(options)) + " enantiomers when chiral centers with " +
        "no specified chirality are systematically varied."
    )

    num_to_keep_initially = thoroughness * max_variants_per_compound
    options = Utils.random_sample(
        options, num_to_keep_initially, ""
    )

    for option in options:
        a_rd_mol = copy.copy(mol.rdkit_mol)
        for idx, chiral in zip(unasignd, option):
            if chiral == "R":
                a_rd_mol.GetAtomWithIdx(idx).SetChiralTag(
                    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW
                )
            elif chiral == "S":
                a_rd_mol.GetAtomWithIdx(idx).SetChiralTag(
                    Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW
                )

        new_mol = MyMol.MyMol(a_rd_mol)

        if not new_mol.remove_bizarre_substruc():
            new_mol.contnr_idx = mol.contnr_idx
            new_mol.name = mol.name
            new_mol.genealogy = mol.genealogy[:]
            new_mol.genealogy.append(
                new_mol.smiles(True) + " (chirality)"
            )
            results.append(new_mol) #, smi
    return results


def enumerate_chiral_molecules(contnrs, max_variants_per_compound, thoroughness, num_processors, multithread_mode, parallelizer_obj):
    """
    Enumerates all possible enantiomers of a molecule. If the chiral of an
    atom is given, that chiral is not varied. Only the chiral of
    unspecified chiral centers is varied.
    """

    if max_variants_per_compound == 0:
        return

    Utils.log("Enumerating all possible enantiomers for all molecules...")

    params = []
    for contnr in contnrs:
        for mol in contnr.mols:
            params.append(tuple([mol, thoroughness, max_variants_per_compound]))
    params = tuple(params)
    tmp = parallelizer_obj.run(params, GetChiral, num_processors, multithread_mode)

    clean = Parallelizer.strip_none(tmp)
    flat = Parallelizer.flatten_list(clean)
    contnr_idxs_prot_failed = Utils.fix_no_prot_generated(contnrs, flat)

    for miss_indx in contnr_idxs_prot_failed:
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
