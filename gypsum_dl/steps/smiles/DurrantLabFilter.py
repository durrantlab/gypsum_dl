"""
This module removes molecules with prohibited substructures, per Durrant-lab
filters.
"""

import __future__

import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl import chem_utils, utils

try:
    from rdkit import Chem
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")

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

# Get the substrings you won't permit (per substring matching)
prohibited_smi_substrs_for_substr = [
    # Let's eliminate ones with common metals too (not really druglike)
    # "[#13]",  # Al
    # "[#23]",  # V
    # "[#26]",  # Fe
    # "[#27]",  # Co
    # "[#29]",  # Cu
    # "[#30]",  # Zn
    # "[#42]",  # Mo
    # "[#48]",  # Cd
    # "[#79]",  # Au
    # "[#82]"   # Pb
    # "[#83]",  # Bi
    "[Al",  # Al
    "[V",  # V
    "[Fe",  # Fe
    "[Co",  # Co
    "[Cu",  # Cu
    "[Zn",  # Zn
    "[Mo",  # Mo
    "[Cd",  # Cd
    "[Au",  # Au
    "[Pb",  # Pb
    "[Bi",  # Bi
]


def durrant_lab_contains_bad_substr(smiles):
    """Determines if a smiles string contains a prohibitive substring. Faster
    than substructure matching.

    :param smiles: The SMILES string.
    :type smiles: A string.
    :return: True if it contains the substring. False otherwise.
    :rtype: boolean
    """

    return any(s in smiles for s in prohibited_smi_substrs_for_substr)


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

    prohibited_substructs = [
        Chem.MolFromSmarts(s) for s in prohibited_smi_substrs_for_substruc
    ]

    # Get the parameters to pass to the parallelizer object.
    params = [[c, prohibited_substructs] for c in contnrs]

    # Run the tautomizer through the parallel object.
    tmp = []
    if parallelizer_obj is None:
        tmp.extend(
            parallel_durrant_lab_filter(c, prohibited_substructs)
            for c, prohibited_substructs in params
        )
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

    # Using this function just to make the changes. Doesn't do energy
    # minimization or anything (as it does later) because max variants
    # and thoroughness maxed out.
    chem_utils.bst_for_each_contnr_no_opt(
        contnrs,
        mols,
        1000,
        1000,  # max_variants_per_compound, thoroughness
    )


def parallel_durrant_lab_filter(contnr, prohibited_substructs):
    """A parallelizable helper function that checks that tautomers do not
       break any nonaromatic rings present in the original object.

    :param contnr: The molecule container.
    :type contnr: MolContainer.MolContainer
    :param prohibited_substructs: A list of the prohibited substructures.
    :type prohibited_substructs: list
    :return: Either the container with bad molecules removed, or a None
      object.
    :rtype: MolContainer.MolContainer | None
    """

    # Replace any molecules that have prohibited substructure with None.
    for mi, m in enumerate(contnr.mols):
        for pattrn in prohibited_substructs:
            if durrant_lab_contains_bad_substr(
                m.orig_smi_deslt
            ) or m.rdkit_mol.HasSubstructMatch(pattrn):
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

                # continue # JDD: this was wrong, wasn't it?
                break  # On to next mol in mols.

    # Now go back and remove those Nones
    contnr.mols = Parallelizer.strip_none(contnr.mols)

    # If there are no molecules, mark this container for deletion.
    return None if len(contnr.mols) == 0 else contnr
