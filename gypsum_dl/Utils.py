# Copyright 2018 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Some helpful utility definitions used throughout the code.
"""

import __future__

import subprocess
import textwrap
import random
import string

def group_mols_by_container_index(mol_lst, job_manager):
    """Take a list of MyMol.MyMol objects, and place them in lists according to
    their associated contnr_idx values. These lists are accessed via
    a dictionary, where they keys are the contnr_idx values
    themselves.

    :param mol_lst: The list of MyMol.MyMol objects.
    :type mol_lst: list
    :param job_manager: The name of the job manager being used.
    :type job_manager: string
    :return: A dictionary, where keys are contnr_idx values and values are
       lists of MyMol.MyMol objects
    :rtype: dict
    """

    # Make the dictionary.
    grouped_results = {}
    for mol in mol_lst:
        if mol is None:
            # Ignore molecules that are None.
            continue

        idx = mol.contnr_idx

        if not idx in grouped_results:
            grouped_results[idx] = []
        grouped_results[idx].append(mol)

    # Remove redundant entries.
    for key in list(grouped_results.keys()):
        grouped_results[key] = list(set(grouped_results[key]))

    # When running in MPI mode, there should be only one container, and it's
    # index won't necessarily match. Let's deal with that here.
    print(">>>>>>>")
    grouped_results = fix_contnr_grpings_when_mpi(grouped_results, mol_lst, job_manager)
    print("<<<<<<<")

    return grouped_results

def random_sample(lst, num, msg_if_cut=""):
    """Randomly selects elements from a list.

    :param lst: The list of elements.
    :type lst: list
    :param num: The number to randomly select.
    :type num: int
    :param msg_if_cut: The message to display if some elements must be ignored
       to construct the list. Defaults to "".
    :param msg_if_cut: str, optional
    :return: A list that contains at most num elements.
    :rtype: list
    """

    try:
        # Remove redundancies.
        lst = list(set(lst))
    except:
        # Because someitems lst element may be unhashable.
        pass

    # Shuffle the list.
    random.shuffle(lst)
    if num < len(lst):
        # Keep the top ones.
        lst = lst[:num]
        if msg_if_cut != "":
            log(msg_if_cut)
    return lst

def log(txt):
    """Prints a message to the screen.

    :param txt: The message to print.
    :type txt: str
    """

    whitespace_before = txt[:len(txt) - len(txt.lstrip())].replace("\t", "    ")
    print(textwrap.fill(
        txt.strip(), width = 80, initial_indent = whitespace_before,
        subsequent_indent = whitespace_before + "    "
    ))

def fix_contnr_grpings_when_mpi(dict_to_fix, mymol_lst, job_manager):
    """In mpi mode, there are often mismatches between dictionaries that map
    container indexes to lists of SMILES or MyMol.MyMol objects, and the actual
    container indexes. This is because of how the containers are divided for
    processing on multiple mpi processors. This function detects this problem
    and fixes it.

    :param dict_to_fix: The dictionary that maps container indexes (right or
       wrong) to SMILES or MyMol.MyMol objects.
    :type dict_to_fix: dict
    :param mymol_lst: A list of the existing containers (with correct
       contnr_idx values).
    :type mymol_lst: list
    :param job_manager: The name of the job manager being used.
    :type job_manager: string
    :return: A dictionary, mapping correct container indexes to the same lists
       (SMILES or MyMol.MyMol).
    :rtype: dict
    """

    return dict_to_fix

    print("dict_to_fix: ", dict_to_fix)
    print("contnr_idxs: ", [m.contnr_idx for m in mymol_lst])

    # When running in MPI mode, there should be only one container, and it's
    # index won't necessarily match. Let's deal with that here.
    first_dict_to_fix_key = dict_to_fix.keys()[0]
    first_contnr_idx = mymol_lst[0].contnr_idx
    if (job_manager == "mpi" and len(dict_to_fix.keys()) == 1 and
        len(set([m.contnr_idx for m in mymol_lst])) == 1 and
        first_dict_to_fix_key != first_contnr_idx):

        print("HIHO")
        dict_to_fix = {first_contnr_idx: dict_to_fix[first_dict_to_fix_key]}
    print("dict_to_fix2: ", dict_to_fix)

    return dict_to_fix

def fnd_contnrs_not_represntd(contnrs, results, job_manager):
    """Identify containers that have no representative elements in results.
    Something likely failed for the containers with no results.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param results: A list of MyMol.MyMol objects.
    :type results: list
    :param job_manager: The name of the job manager being used.
    :type job_manager: string
    :return: A list of integers, the indecies of the contnrs that have no
       associated elements in the results.
    :rtype: list
    """

    # Find ones that don't have any generated. In the context of ionization,
    # for example, this is because sometimes Dimorphite-DL failes to producce
    # valid smiles. In this case, just use the original smiles. Couldn't find
    # a good solution to work around.

    # Get a dictionary of all the input smiles. Keys are indexes, values are
    # smiles.
    idx_to_smi = {}
    for idx in range(0,len(contnrs)):
        contnr = contnrs[idx]
        if not idx in idx_to_smi:
            idx_to_smi[idx] = contnrs[idx].orig_smi_deslt

    # When running in MPI mode, there should be only one container, and it's
    # index won't necessarily match. Let's deal with that here.
    idx_to_smi = fix_contnr_grpings_when_mpi(idx_to_smi, results, job_manager)

    # Now remove from those any that have associated ionized smiles strings.
    # These are represented, and so you don't want to include them in the
    # return.
    for m in results:
        if m.contnr_idx in idx_to_smi:
            del idx_to_smi[m.contnr_idx]

    # Return just the container indexes (the keys).
    return list(idx_to_smi.keys())

def print_current_smiles(contnrs):
    """Prints the smiles of the current containers. Helpful for debugging.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    """

    # For debugging.
    log("    Contents of MolContainers")
    for i, mol_cont in enumerate(contnrs):
        log("\t\tMolContainer #" + str(i) + " (" + mol_cont.name + ")")
        for i, s in enumerate(mol_cont.all_can_noh_smiles()):
            log("\t\t\tMol #" + str(i) + ": " + s)

def exception(msg):
    """Prints an error to the screen and raises an exception.

    :param msg: The error message.
    :type msg: str
    :raises Exception: The error.
    """

    log(msg)
    log("\n" + "=" * 79)
    log("For help with usage:")
    log("\tpython run_gypsum_dl.py --help")
    log("=" * 79)
    log("")
    raise Exception(msg)

def slug(strng):
    """Converts a string to one that is appropriate for a filename.

    :param strng: The input string.
    :type strng: str
    :return: The filename appropriate string.
    :rtype: str
    """

    # See
    # https://stackoverflow.com/questions/295135/turn-a-string-into-a-valid-filename
    valid_chars = "-_.%s%s" % (string.ascii_letters, string.digits)
    if strng == "":
        return "untitled"
    else:
        return "".join([c if c in valid_chars else "_" for c in strng])
