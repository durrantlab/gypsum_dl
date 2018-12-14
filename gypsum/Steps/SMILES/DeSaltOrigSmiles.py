"""Desalts the input SMILES strings. If an input SMILES string contains to
molecule, keep the larger one."""

import __future__

import gypsum.Parallelizer as Parallelizer
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

def desalt_orig_smi(contnrs, num_processors, multithread_mode, parallelizer_obj):
    """If an input molecule has multiple unconnected fragments, this removes
       all but the largest fragment.

    :param contnrs: The molecule containers.
    :type contnrs: list
    :param num_processors: Number of processors to use.
    :type num_processors: int
    :param multithread_mode: The multithread mode.
    :type multithread_mode: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    Utils.log(
        "Desalting all molecules (i.e., keeping only largest fragment)."
    )

    # Desalt each of the molecule containers.
    # @@@@@@@@@@@ JAKE FIX LATER TO MULTI BELOW. JDD: Needs attention? Maybe
    # not. This is probably fast.
    # tmp = mp.MultiThreading(contnrs, num_processors, desalter)
    tmp = [desalter(x) for x in contnrs]

    # Go through each contnr and update the orig_smi_deslt. If we update it,
    # also add a note in the genealogy record.
    tmp = Parallelizer.strip_none(tmp)
    for desalt_mol in tmp:
        idx = desalt_mol.contnr_idx
        cont = contnrs[idx]
        if contnrs[desalt_mol.contnr_idx].orig_smi != desalt_mol.orig_smi:
            desalt_mol.genealogy.append(desalt_mol.orig_smi_deslt + " (desalted)")
            cont.update_orig_smi(desalt_mol.orig_smi_deslt)
        cont.add_mol(desalt_mol)

def desalter(contnr):
    """Desalts molecules in a molecule container.

    :param contnr: The molecule container.
    :type contnr: MolContainer.MolContainer
    :return: A molecule object.
    :rtype: MyMol.MyMol
    """

    # Split it into fragments
    frags = contnr.get_frags_of_orig_smi()

    if len(frags) == 1:
        # It's only got one fragment, so default assumption that
        # orig_smi = orig_smi_deslt is correct.
        return contnr.mol_orig_smi
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
        new_mol.genealogy = contnr.mol_orig_smi.genealogy
        new_mol.makeMolFromSmiles() # Need to update the mol.
        return new_mol
