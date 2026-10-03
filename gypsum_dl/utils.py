"""Some helpful utility definitions used throughout the code."""

import contextlib
import math
import random
import string
import textwrap
from typing import TYPE_CHECKING

from loguru import logger

if TYPE_CHECKING:
    # Annotations only. Importing MolContainer at run time would make utils the
    # head of a utils -> MolContainer -> MyMol -> MolObjectHandling -> utils
    # cycle, which MolObjectHandling had to work around with a function-local
    # import. With this module importing nothing from the package, no such
    # workaround is needed anywhere.
    from gypsum_dl.MolContainer import MolContainer
    from gypsum_dl.MyMol import MyMol


def group_mols_by_container_index(mol_lst):
    """Take a list of MyMol.MyMol objects, and place them in lists according to
    their associated contnr_idx values. These lists are accessed via
    a dictionary, where they keys are the contnr_idx values
    themselves.

    Args:
        mol_lst: The list of MyMol.MyMol objects.

    Returns:
        A dictionary, where keys are `contnr_idx` values and values are lists of
            MyMol.MyMol objects.
    """

    # Make the dictionary.
    grouped_results = {}
    for mol in mol_lst:
        if mol is None:
            # Ignore molecules that are None.
            continue

        idx = mol.contnr_idx

        if idx not in grouped_results:
            grouped_results[idx] = []
        grouped_results[idx].append(mol)

    # Remove redundant entries. dict.fromkeys dedups in first-seen order;
    # set() order depends on PYTHONHASHSEED (MyMol hashes its canonical SMILES
    # string, and CPython randomizes str hashing per invocation), and these
    # lists feed the seeded sampling in random_sample.
    for key in list(grouped_results.keys()):
        grouped_results[key] = list(dict.fromkeys(grouped_results[key]))

    return grouped_results


def random_sample(lst: list, num: int, msg_if_cut: str = ""):
    """Randomly selects elements from a list.

    Args:
        lst: The list of elements.
        num: The number to randomly select.
        msg_if_cut: The message to display if some elements must be ignored
            to construct the list. Defaults to `""`.

    Returns:
        A list that contains at most num elements.
    """

    # Copy before shuffling. The dedup below rebinds lst to a new list, but
    # only when the elements are hashable; when it raises, lst is still the
    # caller's list and random.shuffle would reorder it in place.
    lst = list(lst)

    with contextlib.suppress(TypeError):
        # Remove redundancies. Suppress because sometimes an lst element may
        # be unhashable. dict.fromkeys dedups in first-seen order, whereas
        # set() order depends on PYTHONHASHSEED, which would leave the shuffle
        # below varying between runs that share a random_seed.
        lst = list(dict.fromkeys(lst))

    # Shuffle the list.
    random.shuffle(lst)
    if num < len(lst):
        # Keep the top ones.
        lst = lst[:num]
        if msg_if_cut != "":
            log(msg_if_cut)
    return lst


def log(txt: str, trailing_whitespace: str = "") -> None:
    """Prints a message to the screen and passes it to loguru.

    Args:
        txt: The message to print.
        trailing_whitespace: White space to add to the end of the
            message, after the trim. "" by default.
    """

    # Wrap each line independently so that embedded newlines (e.g. a list of
    # failed SMILES joined with "\n") are preserved instead of being collapsed
    # into a single reflowed paragraph. Long tokens are never split: most
    # messages here embed a SMILES string, and a SMILES broken across lines
    # (or at a hyphen) can no longer be copied back into a tool or grepped
    # for, which is exactly what a user needs when a compound disappears from
    # the output.
    wrapped_lines = []
    for line in txt.split("\n"):
        whitespace_before = line[: len(line) - len(line.lstrip())].replace("\t", "    ")
        wrapped_lines.append(
            textwrap.fill(
                line.strip(),
                width=80,
                initial_indent=whitespace_before,
                subsequent_indent=f"{whitespace_before}    ",
                break_long_words=False,
                break_on_hyphens=False,
            )
        )
    text = "\n".join(wrapped_lines)
    print(text + trailing_whitespace)

    # Also hand the message to loguru, which drops it unless the GYPSUM_DL_LOG
    # environment variables (or enable_logging) turned logging on. No depth
    # offset: loguru filters on the calling module's name, and only messages
    # attributed to a gypsum_dl module are covered by logger.disable, so a
    # caller outside the package would otherwise reach loguru's default stderr
    # sink.
    logger.info(text)


def fnd_contnrs_not_represntd(contnrs: list["MolContainer"], results: list) -> list:
    """Identify containers that have no representative elements in results.
    Something likely failed for the containers with no results.

    Args:
        contnrs: A list of containers (MolContainer.MolContainer).
        results: A list of MyMol.MyMol objects.

    Returns:
        A list of integers, the indecies of the contnrs that have no
            associated elements in the results.
    """

    # Find ones that don't have any generated. In the context of ionization, for
    # example, this is because sometimes Dimorphite-DL failes to producce valid
    # smiles. In this case, just use the original smiles. Couldn't find a good
    # solution to work around.

    # Collect the contnr_idx of every container that produced a result. Key by
    # contnr_idx (not list position): if contnrs was ever filtered, position no
    # longer equals contnr_idx, and results carry contnr_idx.
    represented = {m.contnr_idx for m in results if m is not None}

    # Return the contnr_idx of containers with no representative results.
    return [
        contnr.contnr_idx for contnr in contnrs if contnr.contnr_idx not in represented
    ]


def contnrs_by_idx(
    contnrs: list["MolContainer"],
) -> dict[int, "MolContainer"]:
    """Index a container list by contnr_idx instead of by list position.

    Molecules carry contnr_idx, and so do the failure lists built from them
    (fnd_contnrs_not_represntd). That value equals a container's position in
    contnrs only as long as nothing filters, reorders, or renumbers the list,
    which is not something the steps can see: mpi mode already renumbers every
    container to zero, and it is only safe there because each job holds a
    single container. Looking containers up through this map keeps a step
    correct either way, and turns a mismatch into a KeyError at the lookup
    rather than a wrong container (or an IndexError) much later.

    Args:
        contnrs: A list of containers (MolContainer.MolContainer).

    Returns:
        A dictionary mapping each container's contnr_idx to that container.

    Raises:
        Exception: If two containers share a contnr_idx.
    """

    by_idx: dict[int, "MolContainer"] = {}
    for contnr in contnrs:
        if contnr.contnr_idx in by_idx:
            # Two containers sharing an index makes every regrouping step
            # ambiguous, and a dict would quietly keep whichever came last.
            # Positional indexing hid the same problem by quietly keeping
            # whichever came first.
            exception(
                "Two molecule containers share contnr_idx "
                f"{contnr.contnr_idx} ({by_idx[contnr.contnr_idx].name} and "
                f"{contnr.name}). Container indices must be unique."
            )
        by_idx[contnr.contnr_idx] = contnr

    return by_idx


def carry_over_unrepresented(
    contnrs: list["MolContainer"],
    results: list["MyMol"],
    variant_desc: str,
    genealogy_note: str,
) -> None:
    """Return a container's existing variants to results when a step made none.

    A compound that an enumeration step cannot process should reach the next
    step as it arrived, with a log line and a genealogy entry saying so. Each
    step used to carry its own copy of this fallback, and the copies drifted:
    one step had none, so its failures went unrecorded in the genealogy. The
    carried-over molecules are appended to results rather than written to the
    container because bst_for_each_contnr_no_opt, which every caller runs next,
    repopulates each container from that list.

    Args:
        contnrs: A list of containers (MolContainer.MolContainer).
        results: The step's flat list of new variants. Extended in place.
        variant_desc: What the step produces (e.g., "tautomers"), for the log.
        genealogy_note: The entry appended to each carried-over molecule's
            genealogy.
    """

    contnr_by_idx = contnrs_by_idx(contnrs)
    for miss_indx in fnd_contnrs_not_represntd(contnrs, results):
        failed_contnr = contnr_by_idx[miss_indx]
        log(
            "\tCould not generate valid "
            + variant_desc
            + " for "
            + failed_contnr.orig_smi
            + " ("
            + failed_contnr.name
            + "), so using existing "
            + "(unprocessed) structures."
        )
        for mol in failed_contnr.mols:
            mol.genealogy.append(genealogy_note)
            results.append(mol)


def print_current_smiles(contnrs: list["MolContainer"]) -> None:
    """Prints the smiles of the current containers. Helpful for debugging.

    Args:
        contnrs: A list of containers (MolContainer.MolContainer).
    """

    # For debugging.
    log("    Contents of MolContainers")
    for i, mol_cont in enumerate(contnrs):
        log("\t\tMolContainer #" + str(i) + " (" + mol_cont.name + ")")
        for i, s in enumerate(mol_cont.all_can_noh_smiles()):
            # smiles() reports failure as None, which would otherwise make this
            # debug dump the thing that ends the run.
            log("\t\t\tMol #" + str(i) + ": " + str(s))


def exception(msg: str) -> None:
    """Prints an error to the screen and raises an exception.

    Args:
        msg: The error message.
    """

    log(msg)
    log("\n" + "=" * 79)
    log("For help with usage:")
    log("\tgypsum-dl --help")
    log("=" * 79)
    log("")
    raise Exception(msg)


def energy_for_output(energy: float) -> float | None:
    """Convert a conformer energy into a value fit to publish as a property.

    A conformer whose force field could not be set up carries an infinite
    energy, which is what keeps it from being selected over a conformer that
    was actually scored. That sentinel must not reach the output files, where
    it reads as a measurement and breaks consumers that parse the field as a
    number. Returning None hands the decision to set_rdkit_mol_prop, which
    leaves a None-valued property out of the file entirely.

    Args:
        energy: The conformer energy, in kcal/mol.

    Returns:
        The energy when it is finite, otherwise None.
    """

    return energy if math.isfinite(energy) else None


def describe_energy(energy: float) -> str:
    """Render a conformer energy for the genealogy record.

    Same reasoning as energy_for_output: the genealogy is read by users
    tracking down where a variant came from, so a failed force field should say
    so rather than print an infinity with a unit after it.

    Args:
        energy: The conformer energy, in kcal/mol.

    Returns:
        The energy with its unit, or a note that no energy is available.
    """

    if math.isfinite(energy):
        return f"{energy} kcal/mol"
    return "energy unavailable, force field failed"


def slug(strng: str) -> str:
    """Converts a string to one that is appropriate for a filename.

    Args:
        strng: The input string.

    Returns:
        The filename appropriate string.
    """

    # See
    # https://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename
    if strng == "":
        return "untitled"

    valid_chars = f"-_.{string.ascii_letters}{string.digits}"
    return "".join([c if c in valid_chars else "_" for c in strng])
