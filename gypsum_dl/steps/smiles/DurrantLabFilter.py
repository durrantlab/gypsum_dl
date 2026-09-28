"""
This module removes molecules with prohibited substructures, per Durrant-lab
filters.
"""

import __future__

from functools import cache
from typing import TYPE_CHECKING

import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl import chem_utils, utils

try:
    from rdkit import Chem
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")

if TYPE_CHECKING:
    from gypsum_dl.MolContainer import MolContainer

# Get the substructures you won't permit (per substructure matching, not
# substring matching)
prohibited_smi_substrs_for_substruc = [
    "C=[N-]",
    "[N-]C=[N+]",
    "[nH+]c[n-]",
    "[#7+]~[#7+]",
    "[#7-]~[#7-]",
    "[!#7]~[#7+]~[#7-]~[!#7]",  # Doesn't hit azide.
    # Vina can't process boron anyway...
    "[#5]",  # B
    "O=[PH](=O)([#8])([#8])",  # molvs does odd tautomer: OP(O)(O)=O => O=[PH](=O)(O)O
    "[#7]=C1[#7]=C[#7]C=C1",  # Prevents an odd tautomer sometimes seen with adenine.
    "N=c1cc[#7]c[#7]1",  # Variant of above
    # "[$([NX2H1]),$([NX3H2])]=C[$([OH]),$([O-])]",  # Terminal iminol
    "[$(N)]=C[$([OH]),$([O-])]",  # iminol (including internal)
    "[$(N)]C(=C)[$([OH]),$([O-])]",  # A mistaken amide tautomer that sometimes arises
]

# Metals are not druglike, so every one of them is rejected. The list stops at
# uranium: everything heavier is synthetic and will not appear in a ligand
# file. Metalloids (Si, Ge, As, Sb, Te, Po, At) are deliberately absent, and
# boron has its own [#5] entry in the substructure list above.
metal_element_symbols: list[str] = (
    # Alkali and alkaline earth metals.
    "Li Na K Rb Cs Fr Be Mg Ca Sr Ba Ra "
    # Transition metals.
    "Sc Ti V Cr Mn Fe Co Ni Cu Zn "
    "Y Zr Nb Mo Tc Ru Rh Pd Ag Cd "
    "Hf Ta W Re Os Ir Pt Au Hg "
    # Lanthanides.
    "La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu "
    # Actinides, through uranium.
    "Ac Th Pa U "
    # Post-transition metals.
    "Al Ga In Sn Tl Pb Bi"
).split()

# Get the substrings you won't permit (per substring matching). A metal can
# only be written as a bracketed atom, so keeping the opening bracket is what
# separates indium from iodine and sodium from nitrogen. One symbol overreaches
# slightly: "[K" also catches krypton, which is no more druglike than a metal.
prohibited_smi_substrs_for_substr = [f"[{sym}" for sym in metal_element_symbols]


def durrant_lab_contains_bad_substr(smiles):
    """Determines if a smiles string contains a prohibitive substring. Faster
    than substructure matching.

    :param smiles: The SMILES string.
    :type smiles: A string.
    :return: True if it contains the substring. False otherwise.
    :rtype: boolean
    """

    return any(s in smiles for s in prohibited_smi_substrs_for_substr)


@cache
def get_prohibited_substructs() -> tuple["Chem.Mol", ...]:
    """Compile the prohibited substructure queries once per process.

    The queries used to be built in the parent and handed to the parallelizer
    alongside each container, which pushed them through the task queue (or the
    mpi scatter) once per container and made the filtering depend on RDKit
    preserving SMARTS query features across a pickle round trip. Building them
    inside the worker keeps the serial and parallel paths identical by
    construction, whatever the installed RDKit does.

    Returns:
        The compiled queries, in the order the patterns are declared.
    """

    return tuple(Chem.MolFromSmarts(s) for s in prohibited_smi_substrs_for_substruc)


def durrant_lab_filters(contnrs, num_procs, job_manager, parallelizer_obj):
    """Removes any molecules that contain prohibited substructures, per the
    durrant-lab filters.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: A list.
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param job_manager: The multithred mode to use.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    utils.log("Applying Durrant-lab filters to all molecules...")

    # Get the parameters to pass to the parallelizer object. The queries
    # themselves are compiled in the worker, not sent through the queue.
    params = [[c] for c in contnrs]

    # Run the tautomizer through the parallel object.
    tmp = []
    if parallelizer_obj is None:
        tmp.extend(Parallelizer.run_one(parallel_durrant_lab_filter, i) for i in params)
    else:
        tmp = parallelizer_obj.run(
            params, parallel_durrant_lab_filter, num_procs, job_manager
        )
    # Note that results is a list of containers.

    # Stripping out None values (failed).
    results = Parallelizer.strip_none(tmp)

    # You need to get the molecules as a flat array so you can run it through
    # bst_for_each_contnr_no_opt
    mols = []
    for contnr in results:
        mols.extend(contnr.mols)

    # Also clear contnrs, because they will be re-added using
    # bst_for_each_contnr_no_opt below.
    for contnr in contnrs:
        contnr.mols = []

    # contnrs = results

    # Using this function just to make the changes. It must not drop variants
    # or generate conformers here, which means the cap has to be at least as
    # large as any one container's share of the candidates; no container can
    # hold more than the whole list, so the list length is a cap that keeps
    # pick_lowest_enrgy_mols on its return-everything path whatever
    # max_variants_per_compound the user asked for. A fixed cap silently
    # pruned (and embedded throwaway conformers for) any container that had
    # been allowed to grow past it.
    variant_cap = max(1, len(mols))

    # Carry-over must be off: contnr.mols was just emptied above, so there are
    # no originals left to fall back on. Leaving it on makes a fully filtered
    # container log that its original conformers were kept, contradicting both
    # the per-variant "discarding it" message and gypsum_dl_failed.smi.
    chem_utils.bst_for_each_contnr_no_opt(
        contnrs,
        mols,
        variant_cap,
        variant_cap,  # max_variants_per_compound, thoroughness
        crry_ovr_frm_lst_step_if_no_fnd=False,
    )


def parallel_durrant_lab_filter(contnr: "MolContainer") -> "MolContainer | None":
    """A parallelizable helper function that checks that tautomers do not
       break any nonaromatic rings present in the original object.

    :param contnr: The molecule container.
    :type contnr: MolContainer.MolContainer
    :return: Either the container with bad molecules removed, or a None
      object.
    :rtype: MolContainer.MolContainer | None
    """

    prohibited_substructs = get_prohibited_substructs()

    # Replace any molecules that have prohibited substructure with None.
    for mi, m in enumerate(contnr.mols):
        # A variant whose RDKit mol failed to build can't be substructure
        # matched, and its canonical SMILES can't be generated either, so it
        # gets its own message rather than the one below.
        if m.rdkit_mol is None:
            utils.log(
                "\tA variant generated from "
                + contnr.orig_smi
                + " ("
                + m.name
                + ") has no valid RDKit molecule, so I'm discarding it."
            )
            contnr.mols[mi] = None
            continue

        # The substring test looks at the molecule alone, so it belongs outside
        # the pattern loop; inside, it was re-evaluated once per query. Both
        # tests still short circuit, so a molecule is matched against no more
        # patterns than before.
        prohibited = durrant_lab_contains_bad_substr(m.orig_smi_deslt) or any(
            m.rdkit_mol.HasSubstructMatch(pattrn) for pattrn in prohibited_substructs
        )

        if prohibited:
            utils.log(
                "\t"
                + m.smiles(True)
                + ", a variant generated "
                + "from "
                + contnr.orig_smi
                + " ("
                + m.name
                + "), contains a prohibited substructure, so I'm "
                + "discarding it."
            )

            contnr.mols[mi] = None

    # Now go back and remove those Nones
    contnr.mols = Parallelizer.strip_none(contnr.mols)

    # If there are no molecules, mark this container for deletion.
    return None if len(contnr.mols) == 0 else contnr
