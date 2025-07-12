"""
Desalts the input SMILES strings. If an input SMILES string contains to
molecule, keep the larger one.
"""

from loguru import logger

import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl.models import Molecule


def desalt_orig_smi(
    contnrs, num_procs, job_manager, parallelizer_obj, durrant_lab_filters=False
):
    """If an input molecule has multiple unconnected fragments, this removes
       all but the largest fragment.

    :param contnrs: A list of containers (container.MoleculeContainer).
    :type contnrs: list
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param job_manager: The multiprocess mode.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    logger.info("Desalting all molecules (i.e., keeping only largest fragment).")

    # Desalt each of the molecule containers. This step is very fast, so let's
    # just run it on a single processor always.
    tmp = [desalter(x) for x in contnrs]

    # Go through each contnr and update the orig_smi_deslt. If we update it,
    # also add a note in the genealogy record.
    tmp = Parallelizer.strip_none(tmp)
    for idx in range(len(tmp)):
        desalt_mol = tmp[idx]
        # idx = desalt_mol.contnr_idx
        cont = contnrs[idx]

        if cont.orig_smi != desalt_mol.orig_smi:
            desalt_mol.genealogy.append(f"{desalt_mol.orig_smi_deslt} (desalted)")
            cont.update_orig_smi(desalt_mol.orig_smi_deslt)

        cont.add_mol(desalt_mol)


def desalter(contnr):
    """Desalts molecules in a molecule container.

    :param contnr: The molecule container.
    :type contnr: container.MoleculeContainer
    :return: A molecule object.
    :rtype: Molecule
    """

    # Split it into fragments
    frags = contnr.initial_molecule.fragments

    if len(frags) == 1:
        # It's only got one fragment, so default assumption that
        # orig_smi = orig_smi_deslt is correct.
        return contnr.mol_orig_frm_inp_smi
    logger.debug(
        "Multiple fragments found in "
        + contnr.initial_molecule.canonical_smiles()
        + " ("
        + contnr.initial_molecule.name
        + ")"
    )

    # Find the biggest fragment
    num_heavy_atoms = []
    num_heavy_atoms_to_frag = {}

    for f in frags:
        num = f.GetNumHeavyAtoms()
        num_heavy_atoms.append(num)
        num_heavy_atoms_to_frag[num] = f

    biggest_frag = num_heavy_atoms_to_frag[max(num_heavy_atoms)]

    # Return info about that biggest fragment.
    new_mol = Molecule(biggest_frag)
    new_mol.contnr_idx = contnr.container_id
    new_mol.name = contnr.name
    new_mol.genealogy = contnr.mol.genealogy
    new_mol.make_mol_frm_smiles_sanitze()  # Need to update the mol.
    return new_mol
