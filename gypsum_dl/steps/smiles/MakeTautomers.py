"""This module makes alternate tautomeric states, using MolVS."""

import __future__

import random
from typing import TYPE_CHECKING, TypedDict

import gypsum_dl.MolObjectHandling as MOH
import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl import MyMol, chem_utils, utils

if TYPE_CHECKING:
    # Annotations only, for the two facts builders below.
    from gypsum_dl.MolContainer import MolContainer

try:
    from rdkit import Chem
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")

try:
    from molvs import tautomer
except Exception:
    utils.exception("You need to install molvs and its dependencies.")


class DoubleBondStereo(TypedDict):
    """One double bond whose geometry a tautomer's parent variant specified.

    Recorded as cis or trans relative to a fixed pair of neighbor atoms rather
    than as E or Z, because tautomerization can change the CIP ranks that E
    and Z are defined by while the geometry stays the same.
    """

    begin: int
    end: int
    stereo_atoms: tuple[int, int]
    cis: bool


class ContnrRingFacts(TypedDict):
    """What the aromaticity filter needs to know about the input compound.

    The filter runs once per candidate tautomer, and a compound produces many.
    Shipping the MolContainer to each of those jobs meant pickling every
    variant of the compound once per tautomer of it, when the comparison needs
    one count and one SMILES for the rejection message.
    """

    orig_smi: str
    num_nonaro_rngs: int


class ContnrChiralityFacts(TypedDict):
    """What the chirality filter needs to know about the input compound.

    Separate from ContnrRingFacts for the same reason the two filters are
    separate functions: each job should carry what its own comparison reads and
    nothing else.
    """

    orig_smi: str
    num_specif_chiral_cntrs: int
    num_unspecif_chiral_cntrs: int


def ring_facts(contnr: "MolContainer") -> ContnrRingFacts:
    """Extract the aromaticity filter's reference counts from a container.

    Args:
        contnr: The container describing the input compound.

    Returns:
        The reference facts for that compound.
    """

    return {
        "orig_smi": contnr.orig_smi,
        "num_nonaro_rngs": contnr.num_nonaro_rngs,
    }


def chirality_facts(contnr: "MolContainer") -> ContnrChiralityFacts:
    """Extract the chirality filter's reference counts from a container.

    Args:
        contnr: The container describing the input compound.

    Returns:
        The reference facts for that compound.
    """

    return {
        "orig_smi": contnr.orig_smi,
        "num_specif_chiral_cntrs": contnr.num_specif_chiral_cntrs,
        "num_unspecif_chiral_cntrs": contnr.num_unspecif_chiral_cntrs,
    }


def make_tauts(
    contnrs,
    max_variants_per_compound,
    thoroughness,
    num_procs,
    job_manager,
    let_tautomers_change_chirality,
    parallelizer_obj,
):
    """Generates tautomers of the molecules. Note that some of the generated
    tautomers are not realistic. If you find a certain improbable
    substructure keeps popping up, add it to the list in the
    `prohibited_substructures` definition found with MyMol.py, in the function
    remove_bizarre_substruc().

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: A list.
    :param max_variants_per_compound: To control the combinatorial explosion,
       only this number of variants (molecules) will be advanced to the next
       step.
    :type max_variants_per_compound: int
    :param thoroughness: How many molecules to generate per variant (molecule)
       retained, for evaluation. For example, perhaps you want to advance five
       molecules (max_variants_per_compound = 5). You could just generate five
       and advance them all. Or you could generate ten and advance the best
       five (so thoroughness = 2). Using thoroughness > 1 increases the
       computational expense, but it also increases the chances of finding good
       molecules.
    :type thoroughness: int
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param let_tautomers_change_chirality: Whether to allow tautomers that
      change the number of chiral centers, or the number of those centers
      carrying an assigned stereochemistry.
    :type let_tautomers_change_chirality: bool
    :param job_manager: The multithred mode to use.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    """

    # No need to proceed if there are no max variants.
    if max_variants_per_compound == 0:
        return

    utils.log("Generating tautomers for all molecules...")

    # MolVS treats max_tautomers as a breadth-first depth limiter rather than a
    # cap on the set it hands back: it stops expanding once the set has reached
    # that size, but never truncates. Passing the bare variant cap therefore
    # denies a second expansion pass to any compound whose first pass already
    # filled the quota, putting tautomers that need two proton shifts out of
    # reach. Every other enumeration step budgets generation at thoroughness *
    # max_variants_per_compound and leaves the narrowing to
    # bst_for_each_contnr_no_opt below, so do the same here.
    max_tauts = thoroughness * max_variants_per_compound

    # Create the parameters to feed into the parallelizer object. Pass the
    # molecule itself, as every other enumeration step does, rather than the
    # container plus an index into it: the container holds every variant of the
    # compound, so sending it once per variant made the pickled payload grow
    # with the square of the variant count. The shared props mapping is built
    # once per container.
    reject_chirality_changes = not let_tautomers_change_chirality
    params = []
    for contnr in contnrs:
        props = contnr.contnr_props()
        params.extend(
            (mol, props, max_tauts, reject_chirality_changes) for mol in contnr.mols
        )
    params = tuple(params)

    # Run the tautomizer through the parallel object.
    tmp = []
    if parallelizer_obj is None:
        # MultiThreading with one processor is the same dispatch serial mode
        # uses, and it draws a seed per job; calling run_one directly left
        # this path drawing from whatever generator state happened to be in
        # place, so --random_seed did not reach it.
        tmp = Parallelizer.MultiThreading(params, 1, parallel_make_taut)
    else:
        tmp = parallelizer_obj.run(params, parallel_make_taut, num_procs, job_manager)

    # Flatten the resulting list of lists.
    none_data = tmp
    taut_data = Parallelizer.flatten_list(none_data)

    # Remove bad tautomers. The chirality filter already ran inside
    # parallel_make_taut, the only place the parent variant is still at hand.
    taut_data = tauts_no_break_arom_rngs(
        contnrs, taut_data, num_procs, job_manager, parallelizer_obj
    )

    # Containers left with no tautomers (the worker failed, or the filters
    # above rejected every form) keep their existing structures. This runs
    # after the filters so that it sees what actually survived.
    utils.carry_over_unrepresented(
        contnrs, taut_data, "tautomers", "(WARNING: Unable to generate tautomers)"
    )

    # Keep only the top few compound variants in each container, to prevent a
    # combinatorial explosion.
    chem_utils.bst_for_each_contnr_no_opt(
        contnrs,
        taut_data,
        max_variants_per_compound,
        thoroughness,
        variant_desc="tautomers",
    )


def specified_double_bonds(mol: "Chem.Mol") -> list[DoubleBondStereo]:
    """Record the double-bond geometry a molecule specifies, by atom index.

    MolVS strips the E/Z label from a double bond in every tautomer it returns
    once any one tautomer makes that bond single, the unchanged tautomer
    included. The double-bond step then treats the bond as unspecified and
    enumerates both isomers. Recording the geometry before enumeration lets it
    be put back afterwards.

    Args:
        mol: The exact molecule handed to MolVS, since the record is keyed on
            its atom indices.

    Returns:
        One entry per double bond with a defined geometry.
    """

    ref = Chem.Mol(mol)
    with MOH.legacy_stereo_perception():
        Chem.AssignStereochemistry(ref, cleanIt=True, force=True)

    cis_labels = (Chem.BondStereo.STEREOZ, Chem.BondStereo.STEREOCIS)
    trans_labels = (Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOTRANS)
    specified: list[DoubleBondStereo] = []
    for bond in ref.GetBonds():
        stereo = bond.GetStereo()
        stereo_atoms = tuple(bond.GetStereoAtoms())
        if bond.GetBondType() != Chem.BondType.DOUBLE or len(stereo_atoms) != 2:
            continue
        if stereo not in cis_labels and stereo not in trans_labels:
            continue
        specified.append(
            {
                "begin": bond.GetBeginAtomIdx(),
                "end": bond.GetEndAtomIdx(),
                "stereo_atoms": (stereo_atoms[0], stereo_atoms[1]),
                "cis": stereo in cis_labels,
            }
        )
    return specified


def restore_double_bond_stereo(
    taut: "Chem.Mol", specified: list[DoubleBondStereo]
) -> None:
    """Put the parent's double-bond geometry back on a MolVS tautomer.

    MolVS transforms edit copies of the input, so atom indices carry over.
    A recorded bond is restored only where the tautomer still has a double
    bond between the same atoms, with the same neighbors; where
    tautomerization made it single, it correctly stays unlabeled.

    Args:
        taut: A tautomer returned by MolVS. Modified in place.
        specified: The parent's record, from specified_double_bonds.
    """

    if not specified:
        return

    # Single-bond directions are what the SMILES writer turns into "/" and
    # "\", and MolVS leaves stale ones behind. Clear them all and rebuild them
    # from the restored bond labels, so the two cannot disagree.
    for bond in taut.GetBonds():
        bond.SetBondDir(Chem.BondDir.NONE)
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            bond.SetStereo(Chem.BondStereo.STEREONONE)

    for spec in specified:
        bond = taut.GetBondBetweenAtoms(spec["begin"], spec["end"])
        if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
            continue
        stereo_begin, stereo_end = spec["stereo_atoms"]
        if bond.GetBeginAtomIdx() != spec["begin"]:
            stereo_begin, stereo_end = stereo_end, stereo_begin
        if (
            taut.GetBondBetweenAtoms(stereo_begin, bond.GetBeginAtomIdx()) is None
            or taut.GetBondBetweenAtoms(stereo_end, bond.GetEndAtomIdx()) is None
        ):
            continue
        bond.SetStereoAtoms(stereo_begin, stereo_end)
        bond.SetStereo(
            Chem.BondStereo.STEREOCIS if spec["cis"] else Chem.BondStereo.STEREOTRANS
        )

    Chem.SetDoubleBondNeighborDirections(taut)


def parallel_make_taut(
    mol: MyMol.MyMol,
    props: MyMol.ContnrProps,
    max_tauts: int,
    reject_chirality_changes: bool = False,
) -> list[MyMol.MyMol] | None:
    """Makes alternate tautomers for a given molecule. This is the function
       that gets fed into the parallelizer.

    :param mol: The molecule (variant) to tautomerize.
    :type mol: MyMol.MyMol
    :param props: The container-level fields describing the input compound the
       molecule is a variant of.
    :type props: MyMol.ContnrProps
    :param max_tauts: The size at which MolVS stops expanding the tautomer
       set. It can return more than this many forms, since it does not discard
       any it has already built, and the caller trims to
       max_variants_per_compound afterwards regardless.
    :type max_tauts: int
    :param reject_chirality_changes: Whether to drop tautomers whose chiral
       center counts differ from those of mol. See
       parallel_check_chiral_centers.
    :type reject_chirality_changes: bool
    :return: A list of MyMol.MyMol objects, containing the alternate
        tautomeric forms.
    :rtype: list
    """

    # Create a temporary RDKit mol object, since that's what MolVS works with.
    # TODO: There should be a copy function
    m = MyMol.MyMol(mol.smiles()).rdkit_mol

    # Make sure it's not None.
    if m is None:
        utils.log(
            "\tCould not generate tautomers for "
            + props["orig_smi"]
            + ". I'm deleting it."
        )
        return

    # For tautomers to work, you need to not have any explicit hydrogens.
    m = Chem.RemoveHs(m)
    if m is None:
        return None

    # Molecules should be kekulized already, but let's double check that.
    # Because MolVS requires kekulized input.
    try:
        Chem.Kekulize(m)
    except Exception:
        return None
    m = MOH.check_sanitization(m)
    if m is None:
        return None

    # Stop expanding once the set reaches max_tauts. Note that another batch
    # could add more, and MolVS itself can overshoot, so you'll need to trim to
    # max_variants_per_compound later. But this could at least help prevent the
    # combinatorial explosion at this stage.
    specified = specified_double_bonds(m)
    enum = tautomer.TautomerEnumerator(max_tautomers=max_tauts)
    tauts_rdkit_mols = enum.enumerate(m)
    for taut_rdkit_mol in tauts_rdkit_mols:
        restore_double_bond_stereo(taut_rdkit_mol, specified)

    # Make all those tautomers into MyMol objects.
    tauts_mols = [MyMol.MyMol(m) for m in tauts_rdkit_mols]

    # Keep only those that have reasonable substructures.
    tauts_mols = [t for t in tauts_mols if t.remove_bizarre_substruc() == False]

    # If there's more than one, let the user know that.
    if len(tauts_mols) > 1:
        utils.log("\t" + str(mol.smiles(True)) + " has tautomers.")

    # The chirality reference is the variant being tautomerized, not the input
    # compound the container describes. Ionization runs first and can
    # legitimately add or remove a stereocenter (protonating a tertiary amine
    # with three distinct substituents creates one at [NH+]), and charging that
    # change to tautomerization rejected every form of the variant, the
    # unchanged one included.
    parent_facts: ContnrChiralityFacts = {
        "orig_smi": props["orig_smi"],
        "num_specif_chiral_cntrs": len(mol.chiral_cntrs_only_asignd()),
        "num_unspecif_chiral_cntrs": len(mol.chiral_cntrs_w_unasignd()),
    }

    # Now collect the final results.
    results = []

    for tm in tauts_mols:
        tm.inherit_contnr_props(props)
        tm.genealogy = mol.genealogy[:]
        tm.name = mol.name

        if (
            reject_chirality_changes
            and parallel_check_chiral_centers(tm, parent_facts) is None
        ):
            continue

        if tm.smiles() != mol.smiles():
            tm.genealogy.append(f"{tm.smiles(True)} (tautomer)")

        results.append(tm)

    return results


def tauts_no_break_arom_rngs(
    contnrs, taut_data, num_procs, job_manager, parallelizer_obj
):
    """For a given molecule, the number of aromatic rings should never change
       regardless of tautization, ionization, etc. Any taut whose ring
       aromaticity differs from the original is unlikely to be worth pursuing.
       So remove it.

       The test is run on nonaromatic ring counts rather than aromatic ones.
       Tautomerization shifts protons and bond orders but never opens or closes
       a ring, so the total ring count is fixed across the comparison, which
       makes the two counts equivalent tests.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: A list.
    :param taut_data: A list of MyMol.MyMol objects.
    :type taut_data: list
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param job_manager: The multithred mode to use.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    :return: A list of MyMol.MyMol objects, with certain bad ones removed.
    :rtype: list
    """

    # You need to group the taut_data by container to pass it to the
    # paralleizer. Build an index map so a taut that matches no container is
    # skipped rather than paired with a stale (or unbound) container. Each
    # container's facts are extracted once and shared by all of its tautomers,
    # so the job payload no longer carries the container itself.
    facts_by_idx: dict[int, ContnrRingFacts] = {
        contnr_idx: ring_facts(contnr)
        for contnr_idx, contnr in utils.contnrs_by_idx(contnrs).items()
    }
    params = []
    for taut_mol in taut_data:
        facts = facts_by_idx.get(taut_mol.contnr_idx)
        if facts is None:
            continue
        params.append((taut_mol, facts))
    params = tuple(params)

    # Run it through the parallelizer to remove non-aromatic rings.

    tmp = []
    if parallelizer_obj is None:
        tmp = Parallelizer.MultiThreading(params, 1, parallel_check_nonarom_rings)
    else:
        tmp = parallelizer_obj.run(
            params, parallel_check_nonarom_rings, num_procs, job_manager
        )

    # Stripping out None values (failed).
    return Parallelizer.strip_none(tmp)


def tauts_no_elim_chiral(contnrs, taut_data, num_procs, job_manager, parallelizer_obj):
    """Unfortunately, molvs sees removing chiral specifications as being a
       distinct taut. I imagine there are cases where tautization could
       remove a chiral center, but I think these cases are rare. To compensate
       for the error in other folk's code, let's require that isomerization
       leave both the number of chiral centers and the number of them carrying
       an assigned stereochemistry unchanged. The second count is what rejects
       the molvs artifact: dropping a specification leaves the center in place,
       so the total by itself does not move.

       This compares each tautomer against the input compound, which is wrong
       for tautomers of an ionized variant whose stereocenter count ionization
       changed. make_tauts therefore applies the same check inside
       parallel_make_taut, against the parent variant, and does not call this.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param taut_data: A list of MyMol.MyMol objects.
    :type taut_data: list
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param job_manager: The multithred mode to use.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    :return: A list of MyMol.MyMol objects, with certain bad ones removed.
    :rtype: list
    """

    # You need to group the taut_data by contnr to pass to paralleizer. Build
    # an index map so a taut that matches no container is skipped rather than
    # paired with a stale (or unbound) container. Each container's facts are
    # extracted once and shared by all of its tautomers, so the job payload no
    # longer carries the container itself.
    facts_by_idx: dict[int, ContnrChiralityFacts] = {
        contnr_idx: chirality_facts(contnr)
        for contnr_idx, contnr in utils.contnrs_by_idx(contnrs).items()
    }
    params = []
    for taut_mol in taut_data:
        facts = facts_by_idx.get(taut_mol.contnr_idx)
        if facts is None:
            continue
        params.append((taut_mol, facts))
    params = tuple(params)

    # Run it through the parallelizer.
    tmp = []
    if parallelizer_obj is None:
        tmp = Parallelizer.MultiThreading(params, 1, parallel_check_chiral_centers)
    else:
        tmp = parallelizer_obj.run(
            params, parallel_check_chiral_centers, num_procs, job_manager
        )

    # Stripping out None values
    return [x for x in tmp if x != None]


def parallel_check_nonarom_rings(
    taut: MyMol.MyMol, facts: ContnrRingFacts
) -> MyMol.MyMol | None:
    """A parallelizable helper function that checks that tautomers have the
       same ring aromaticity as the original object. The test is symmetric: a
       tautomer that makes a nonaromatic ring aromatic is rejected alongside
       one that dearomatizes an aromatic ring.

    :param taut: The tautomer to evaluate.
    :type taut: MyMol.MyMol
    :param facts: The original compound's reference counts. See
       ContnrRingFacts.
    :type facts: ContnrRingFacts
    :return: Either the tautomer or a None object.
    :rtype: MyMol.MyMol | None
    """

    # How many nonaromatic rings in the original smiles?
    num_nonaro_rngs_orig = facts["num_nonaro_rngs"]

    # Note that a ring counts as nonaromatic here if any one of its atoms is
    # nonaromatic, applied the same way on both sides of the comparison.
    get_idxs_of_nonaro_rng_atms = len(taut.get_idxs_of_nonaro_rng_atms())
    if get_idxs_of_nonaro_rng_atms == num_nonaro_rngs_orig:
        # Same number of nonaromatic rings as original molecule. Saves the
        # good ones.
        return taut
    else:
        utils.log(
            "\t"
            + str(taut.smiles(True))
            + ", a tautomer generated "
            + "from "
            + facts["orig_smi"]
            + " ("
            + taut.name
            + "), changed the number of non-aromatic rings, so I'm discarding it."
        )


def parallel_check_chiral_centers(
    taut: MyMol.MyMol, facts: ContnrChiralityFacts
) -> MyMol.MyMol | None:
    """A parallelizable helper function that checks that tautomers do not break
       any chiral centers in the original molecule.

       Two counts describe a molecule's chirality, and a tautomer has to match
       the original on both. The total number of chiral centers catches a
       transformation that created or destroyed a stereocenter. The number of
       those centers carrying an assignment catches the MolVS artifact this
       filter exists for (see tauts_no_elim_chiral): a form differing from the
       input only in having dropped a chiral specification. That form has the
       same total, so comparing totals alone would admit it, and
       EnumerateChiralMols runs later in prepare_smiles and expands every
       unassigned center into both R and S, which would put the enantiomer the
       input ruled out into the output.

    :param taut: The tautomer to evaluate.
    :type taut: MyMol.MyMol
    :param facts: The original compound's reference counts. See
       ContnrChiralityFacts.
    :type facts: ContnrChiralityFacts
    :return: Either the tautomer or a None object.
    :rtype: MyMol.MyMol | None
    """

    # chiral_cntrs_w_unasignd reports assigned and unassigned centers alike, so
    # its length is the total; num_unspecif_chiral_cntrs holds that length
    # despite its name. Adding it to the assigned count would tally every
    # assigned center twice.
    num_chiral_cntrs_orig = facts["num_unspecif_chiral_cntrs"]
    num_assignd_chiral_cntrs_orig = facts["num_specif_chiral_cntrs"]

    num_chiral_cntrs_taut = len(taut.chiral_cntrs_w_unasignd())
    num_assignd_chiral_cntrs_taut = len(taut.chiral_cntrs_only_asignd())

    if (
        num_chiral_cntrs_taut == num_chiral_cntrs_orig
        and num_assignd_chiral_cntrs_taut == num_assignd_chiral_cntrs_orig
    ):
        # Same chirality as the original molecule. Save this good one.
        return taut

    rejection_prefix = (
        "\t"
        + facts["orig_smi"]
        + " ==> "
        + str(taut.smiles(True))
        + " (tautomer transformation on "
        + taut.name
        + ") "
    )

    if num_chiral_cntrs_taut != num_chiral_cntrs_orig:
        utils.log(
            rejection_prefix
            + "changed the molecules total number of "
            + "chiral centers from "
            + str(num_chiral_cntrs_orig)
            + " to "
            + str(num_chiral_cntrs_taut)
            + ", so I'm deleting it."
        )
    else:
        utils.log(
            rejection_prefix
            + "kept all "
            + str(num_chiral_cntrs_orig)
            + " chiral centers but changed how many of them carry an assigned "
            + "stereochemistry, from "
            + str(num_assignd_chiral_cntrs_orig)
            + " to "
            + str(num_assignd_chiral_cntrs_taut)
            + ", so I'm deleting it."
        )

    return None
