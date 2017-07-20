import subprocess
import textwrap
import random

def group_mols_by_container_index(mol_lst):
    """
    Take a list of MyMol.MyMol objects, and place them in lists according to
    their associated contnr_idx values. These lists are accessed via
    a dictionary, where they keys are the contnr_idx values
    themselves.

    :param [MyMol.MyMol] mol_list: The list of MyMol.MyMol objects.

    :returns: a dictionary, where keys are contnr_idx values and
              values are lists of MyMol.MyMol objects.
    :rtype: :class:`str` ???
    """

    grouped_results = {}
    for mol in mol_lst:
        idx = mol.contnr_idx

        if not idx in grouped_results:
            grouped_results[idx] = []
        grouped_results[idx].append(mol)

    for key in grouped_results.keys():
        grouped_results[key] = list(set(grouped_results[key]))

    return grouped_results

def random_sample(lst, num, msg_if_cut=""):
    """
    Randomly selets elements from a list.

    :param [any] lst: The list of elements.

    :param int num: The number to randomly select.

    :param str msg_if_cut: The message to display if some elements must be
               ignored to construct the list.

    :returns: a list that contains at most num elements.
    :rtype: :class:`str` ???
    """
    
    try:
        lst = list(set(lst))
    except:
        # Because someitems lst element may be unhashable.
        pass
    
    random.shuffle(lst)
    if num < len(lst):
        lst = lst[:num]
        if msg_if_cut != "":
            log(msg_if_cut)
    return lst

def log(txt):
    """
    Prints a message to the screen.

    :param str txt: The message to print.
    """

    whitespace_before = txt[:len(txt) - len(txt.lstrip())].replace("\t", "    ")
    print(textwrap.fill(
        txt.strip(), width = 80, initial_indent = whitespace_before,
        subsequent_indent = whitespace_before + "    "
    ))

def runit(cmd):
    """
    Runs a command and returns the output as a list.

    :param str cmd: The command.

    :returns: A list containing the output, each line as a separate item
                in the list.
    :rtype: :class:`list` ???
    """
    
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, stdin=subprocess.PIPE,
                            shell = True)
    
    out, err = p.communicate()

    return out.decode().split("\n")  # So returns a list of lines.


# See http://stackoverflow.com/questions/952914/making-a-flat-list-out-of-list-of-lists-in-python
def flatten(list_of_lists):
    """
    Flattens a list of lists to just a list.
    
    :param list list_of_lists: The list of lists.
    
    :returns: A list, flattened.
    :rtype: :class:`list` ???
    """

    return [item for sublist in list_of_lists for item in sublist]

def contnrs_no_touchd(self, results):
    """
    Identify contnrs that have no representative elements in results.

    :param ConfGenerator self: The associated ConfGenerator object.

    :param [MyMol.MyMol] results: A list of MyMol.MyMol objects.

    :returns: a list of integers, the indecies of the contnrs that have
              no associated elements in results.
    :rtype: :class:`str` ???
    """
    
    # Find ones that don't have any generated. This is because sometimes
    # obabel failes to producce valid smiles. In this case, just use the
    # original smiles. Couldn't find a good solution to work around.
    
    # Get a dictionary of all the input smiles.
    idx_to_smi = {}
    for contnr in self.contnrs:
        if not contnr.contnr_idx in idx_to_smi:
            idx_to_smi[contnr.contnr_idx] = self.contnrs[contnr.contnr_idx].orig_smi_deslt
    
    # Now remove from those any that have associated protonated smiles strings
    # from obabel.
    for m in results:
        if m.contnr_idx in idx_to_smi:
            del idx_to_smi[m.contnr_idx]
    
    return idx_to_smi.keys()