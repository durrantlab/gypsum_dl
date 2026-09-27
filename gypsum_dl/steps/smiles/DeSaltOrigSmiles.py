"""
Desalts the input SMILES strings. If an input SMILES string contains to
molecule, keep the larger one.
"""

import __future__

import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl import MyMol, chem_utils, utils

try:
    from rdkit import Chem
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")


def desalt_orig_smi(
    contnrs, num_procs, job_manager, parallelizer_obj, durrant_lab_filters=False
):
    """If an input molecule has multiple unconnected fragments, this removes
       all but the largest fragment.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param job_manager: The multiprocess mode.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    utils.log("Desalting all molecules (i.e., keeping only largest fragment).")

    # Desalting is very fast, so always run it on a single processor. Still
    # route each container through run_one: fragmenting and sanitizing can
    # raise (a fragment that will not sanitize apart from its counter-ion, for
    # instance), and desalting is the first step, so an unguarded exception
    # here aborts the run before any output exists. Every other step already
    # drops a single failed molecule this way.
    for cont in contnrs:
        desalt_mol = Parallelizer.run_one(desalter, (cont,))

        if desalt_mol is None or desalt_mol.rdkit_mol is None:
            # Either desalter raised, or the largest fragment did not survive
            # sanitization. Fall back to the molecule as supplied, which did
            # sanitize, rather than letting a molecule with no rdkit_mol reach
            # the later steps and the writers.
            utils.log(
                "\tWARNING: Could not desalt "
                + cont.orig_smi
                + " ("
                + cont.name
                + "). Keeping the molecule as given."
            )
            orig_mol = cont.mol_orig_frm_inp_smi
            if not orig_mol.genealogy:
                orig_mol.genealogy.append(f"{cont.orig_smi} (source)")
            cont.add_mol(orig_mol)
            continue

        # Update the orig_smi_deslt. If we update it, also add a note in the
        # genealogy record.
        if cont.orig_smi != desalt_mol.orig_smi:
            desalt_mol.genealogy.append(f"{desalt_mol.orig_smi_deslt} (desalted)")
            cont.update_orig_smi(desalt_mol.orig_smi_deslt)

            # update_orig_smi rebuilds mol_orig_frm_inp_smi from the desalted
            # SMILES, and the replacement starts with an empty genealogy.
            # Later steps (ionization above all) seed each variant's record
            # from that object, so the record has to survive the rebuild or
            # the provenance of every desalted variant is lost.
            cont.mol_orig_frm_inp_smi.genealogy = desalt_mol.genealogy[:]

        cont.add_mol(desalt_mol)


def desalter(contnr):
    """Desalts molecules in a molecule container.

    :param contnr: The molecule container.
    :type contnr: MolContainer.MolContainer
    :return: A molecule object.
    :rtype: MyMol.MyMol
    """

    # Desalting is the first step every input passes through, and both
    # branches below derive their genealogy from the pristine molecule, so
    # stamp the input SMILES here rather than in each branch. Guarded so that
    # a second call cannot record the source twice.
    orig_mol = contnr.mol_orig_frm_inp_smi
    if not orig_mol.genealogy:
        orig_mol.genealogy.append(f"{contnr.orig_smi} (source)")

    # Split it into fragments
    frags = contnr.get_frags_of_orig_smi()

    if len(frags) == 1:
        # It's only got one fragment, so default assumption that
        # orig_smi = orig_smi_deslt is correct.
        return orig_mol
    utils.log(
        "\tMultiple fragments found in " + contnr.orig_smi + " (" + contnr.name + ")"
    )

    # Find the biggest fragment. max() returns the first maximal element, so a
    # tie between equal-sized fragments is broken by the order they appear in
    # the input SMILES, which makes the choice reproducible and inspectable.
    biggest_frag = max(frags, key=lambda f: f.GetNumHeavyAtoms())

    # Return info about that biggest fragment.
    new_mol = MyMol.MyMol(biggest_frag)
    new_mol.contnr_idx = contnr.contnr_idx
    new_mol.name = contnr.name
    new_mol.genealogy = orig_mol.genealogy[:]
    new_mol.make_mol_frm_smiles_sanitze()  # Need to update the mol.
    return new_mol
