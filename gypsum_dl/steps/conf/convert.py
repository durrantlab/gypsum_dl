"""
A module to so the 2D to 3D conversion, though the actual code for that
conversion is in Molecule.make_first_3d_conf_no_min()
"""

from typing import TYPE_CHECKING

from loguru import logger

import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl import chem_utils

if TYPE_CHECKING:
    from gypsum_dl import Molecule, MoleculeContainer


def convert_2d_to_3d(
    contnrs: list["MoleculeContainer"],
    max_variants_per_compound: int,
    thoroughness: int,
    num_procs: int,
    job_manager: str,
    parallelizer_obj: object,
) -> None:
    """Converts the 1D smiles strings into 3D small-molecule models.

    Args:
        contnrs: A list of containers (container.MoleculeContainer).
        max_variants_per_compound: To control the combinatorial explosion,
            only this number of variants (molecules) will be advanced to the next
            step.
        thoroughness: How many molecules to generate per variant (molecule)
            retained, for evaluation. For example, perhaps you want to advance five
            molecules (max_variants_per_compound = 5). You could just generate five
            and advance them all. Or you could generate ten and advance the best
            five (so thoroughness = 2). Using thoroughness > 1 increases the
            computational expense, but it also increases the chances of finding good
            molecules.
        num_procs: The number of processors to use.
        job_manager: The multithred mode to use.
        parallelizer_obj: The Parallelizer object.
    """

    logger.info("Converting all molecules to 3D structures.")

    # Make the inputs to pass to the parallelizer.
    params = []
    for contnr in contnrs:
        params.extend((mol,) for mol in contnr.mols)
    params = tuple(params)

    # Run the parallelizer
    tmp = []
    if parallelizer_obj is None:
        tmp.extend(parallel_make_3d(i[0]) for i in params)
    else:
        tmp = parallelizer_obj.run(params, parallel_make_3d, num_procs, job_manager)
    # Remove and Nones from the output, which represent failed molecules.
    clear = Parallelizer.strip_none(tmp)

    # Keep only the top few compound variants in each container, to prevent a
    # combinatorial explosion.
    chem_utils.bst_for_each_contnr_no_opt(
        contnrs, clear, max_variants_per_compound, thoroughness, False
    )


def parallel_make_3d(mol: "Molecule") -> "Molecule | None":
    """Does the 2D to 3D conversion. Meant to run within parallelizer.

    Args:
        mol: The molecule to be converted.

    Returns:
        A Molecule object with the 3D coordinates inside, or None if
            it fails.
    """

    # Initially assume you won't show an error message.
    show_error_msg = False

    if mol.rdkit_mol is None:
        # The rdkit mol is None. Something's gone wrong. Show an error
        # message.
        show_error_msg = True
    elif not mol.remove_bizarre_substruc():
        # Check if it has strange substructures.

        # Perform the conversion.
        mol.make_first_3d_conf_no_min()

        # If there are some conformations, make note of that in the
        # genealogy record.
        if len(mol.conformers) > 0:
            mol.genealogy.append(f"{mol.smiles(True)} (3D coordinates assigned)")
            return mol
        else:
            # No conformers? Show an error. Something's gone wrong.
            show_error_msg = True

    if show_error_msg:
        # Something's gone wrong, so show this error.
        logger.warning(
            "Could not generate 3D geometry for "
            + str(mol.smiles())
            + " ("
            + mol.name
            + "). Molecule "
            + "discarded."
        )

    # If you get here, something's gone wrong...
    return None
