"""Some helpful utility definitions used throughout the code."""

import contextlib
import random
import string
import textwrap

from gypsum_dl import MolContainer, MyMol


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
    """Prints a message to the screen.

    Args:
        txt: The message to print.
        trailing_whitespace: White space to add to the end of the
            message, after the trim. "" by default.
    """

    # Wrap each line independently so that embedded newlines (e.g. a list of
    # failed SMILES joined with "\n") are preserved instead of being collapsed
    # into a single reflowed paragraph.
    wrapped_lines = []
    for line in txt.split("\n"):
        whitespace_before = line[: len(line) - len(line.lstrip())].replace(
            "\t", "    "
        )
        wrapped_lines.append(
            textwrap.fill(
                line.strip(),
                width=80,
                initial_indent=whitespace_before,
                subsequent_indent=f"{whitespace_before}    ",
            )
        )
    print("\n".join(wrapped_lines) + trailing_whitespace)


def fnd_contnrs_not_represntd(contnrs: list[MolContainer], results: list) -> list:
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
    return [contnr.contnr_idx for contnr in contnrs if contnr.contnr_idx not in represented]


def print_current_smiles(contnrs: list[MolContainer]) -> None:
    """Prints the smiles of the current containers. Helpful for debugging.

    Args:
        contnrs: A list of containers (MolContainer.MolContainer).
    """

    # For debugging.
    log("    Contents of MolContainers")
    for i, mol_cont in enumerate(contnrs):
        log("\t\tMolContainer #" + str(i) + " (" + mol_cont.name + ")")
        for i, s in enumerate(mol_cont.all_can_noh_smiles()):
            log("\t\t\tMol #" + str(i) + ": " + s)


def exception(msg: str) -> None:
    """Prints an error to the screen and raises an exception.

    Args:
        msg: The error message.
    """

    log(msg)
    log("\n" + "=" * 79)
    log("For help with usage:")
    log("\tpython run_gypsum_dl.py --help")
    log("=" * 79)
    log("")
    raise Exception(msg)


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