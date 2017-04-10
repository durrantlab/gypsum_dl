from ... import Utils
from ... import ChemUtils
from ... import MyMol
import itertools
import copy
from ... import Multiprocess as mp
import sys
import random

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)

def GetDoubleBonded(mol, params):

    # Get all double bonds that don't have defined stereochemistry
    unasignd = mol.get_double_bonds_without_stereochemistry()

    # Throw out any bond that is in a small ring
    unasignd = [i for i in unasignd
                if not mol.GetBondWithIdx(i).IsInRingSize(3)]
    unasignd = [i for i in unasignd
                if not mol.GetBondWithIdx(i).IsInRingSize(4)]
    unasignd = [i for i in unasignd
                if not mol.GetBondWithIdx(i).IsInRingSize(5)]
    unasignd = [i for i in unasignd
                if not mol.GetBondWithIdx(i).IsInRingSize(6)]
    unasignd = [i for i in unasignd
                if not mol.GetBondWithIdx(i).IsInRingSize(7)]

    # Throw out any bond that has an atom that only participate in
    # that one bond (terminal alkene).
    bonds_to_use = []
    for bond_index in unasignd:
        b = mol.GetBondWithIdx(bond_index)

        a1 = b.GetBeginAtom()
        nb1 = a1.GetBonds()
        if len(nb1) == 1:
            # The only bond is the one you already know about. So
            # don't save.
            continue

        a2 = b.GetEndAtom()
        nb2 = a2.GetBonds()
        if len(nb2) == 1:
            # The only bond is the one you already know about. So
            # don't save.
            continue

        bi1 = [b.GetIdx() for b in a1.GetBonds()]
        bi1.remove(bond_index)

        bi2 = [b.GetIdx() for b in a2.GetBonds()]
        bi2.remove(bond_index)

        bonds_to_use.append((bi1[0], bi2[0]))


    # Get all possible chiral assignments.
    num = len(bonds_to_use)
    if num == 0:
        # There are no unspecified chiral centers, so just keep
        # existing.
        return mol
    elif num == 1:
        options = [[True], [False]]
    else:
        starting = [[True], [False]]
        options = [[True], [False]]
        for i in range(num - 1):
            options = list(itertools.product(options, starting))
            options = [list(itertools.chain(c[0], c[1]))
                       for c in options]


    Utils.log(
        "\t" + mol.smiles(True) + " has " + str(len(options)) +
        " double bonds with unspecified stereochemistry."
    )

    options = Utils.random_sample(
        options, params["max_variants_per_compound"], ""
    )

    for option in options:

        a_rd_mol = copy.copy(mol.rdkit_mol)
        for idxs, chiral in zip(bonds_to_use, option):
            i1, i2 = idxs

            a_rd_mol.GetBondWithIdx(i1).SetBondDir(
                Chem.BondDir.ENDUPRIGHT
            )
            if chiral == True:
                a_rd_mol.GetBondWithIdx(i2).SetBondDir(
                    Chem.BondDir.ENDUPRIGHT
                )
            elif chiral == False:
                a_rd_mol.GetBondWithIdx(i2).SetBondDir(
                    Chem.BondDir.ENDDOWNRIGHT
                )

        new_mol = MyMol.MyMol(a_rd_mol)

        if new_mol.can_smi != False:
            # Sometimes you get an error if there's a bad structure
            # otherwise.

            new_mol.AssignStereochemistry()
            new_mol.contnr_idx = mol.contnr_idx


            new_mol.name = mol.name
            new_mol.genealogy = mol.genealogy[:]
            new_mol.genealogy.append(
                new_mol.smiles(True) + " (cis-trans isomerization)"
            )

            return new_mol #, smi))


def enumerate_double_bonds(self):
    """
    Enumerates all possible cis-trans isomers. If the stereochemistry of a
    double bond is specified, it is not varied. All unspecified double bonds
    are varied.
    """

    if self.params["max_variants_per_compound"] == 0:
        return

    Utils.log(
        "Enumerating all possible cis-trans isomers for all molecules..."
    )

    params = []
    for contnr in self.contnrs:
        for mol in contnr.mols:
            params.append((mol, self.params))

    tmp = mp.MultiThreading(
        params, self.params["num_processors"], GetDoubleBonded
    )

    ChemUtils.bst_for_each_contnr_no_opt(self, tmp)
