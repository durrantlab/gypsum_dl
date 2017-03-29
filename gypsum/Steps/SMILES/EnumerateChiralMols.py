from ... import Utils
from ... import ChemUtils
from ... import mp_queue as mp
from ... import MyMol
import copy
import sys
import itertools
import random

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)


def GetChiral(mol, params):

    # Get all chiral centers that aren't assigned.
    unasignd = [p[0] for p in mol.chiral_cntrs_w_unasignd()
                if p[1] == "?"]
    num = len(unasignd)

    # Get all possible chiral assignments. If chiral is specified,
    # retain it.
    if num == 0:
        # There are no unspecified chiral centers, so just keep
        # existing.
        return mol
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

    num_to_keep_initially = params["thoroughness"] * params["max_variants_per_compound"]
    options = Utils.random_sample(
        options, num_to_keep_initially, ""
    )

    for option in options:
        # print mol.smiles(), mol.name, option, "*"
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

        if not new_mol.crzy_substruc():
            new_mol.contnr_idx = mol.contnr_idx
            new_mol.name = mol.name
            new_mol.genealogy = mol.genealogy[:]
            new_mol.genealogy.append(
                new_mol.smiles(True) + " (chirality)"
            )
            return new_mol #, smi


def enumerate_chiral_molecules(self):
    """
    Enumerates all possible enantiomers of a molecule. If the chiral of an
    atom is given, that chiral is not varied. Only the chiral of
    unspecified chiral centers is varied.
    """

    if self.params["max_variants_per_compound"] == 0:
        return

    Utils.log("Enumerating all possible enantiomers for all molecules...")

    params = []
    for contnr in self.contnrs:
        for mol in contnr.mols:
            params.append((mol, self.params))

    tmp = mp.MultiThreading(params, self.params["num_processors"], GetChiral)

    contnr_indx_no_touch = Utils.contnrs_no_touchd(self, tmp)

    for miss_indx in contnr_indx_no_touch:
        Utils.log(
            "\tCould not generate valid enantiomers for " +
            self.contnrs[miss_indx].orig_smi + " (" +
            self.contnrs[miss_indx].name + "), so using existing " +
            "(unprocessed) structures.")
        for mol in self.contnrs[miss_indx].mols:
            mol.genealogy.append("(WARNING: Unable to generate enantiomers)")
            tmp.append(mol)

    ChemUtils.bst_for_each_contnr_no_opt(self, tmp)
