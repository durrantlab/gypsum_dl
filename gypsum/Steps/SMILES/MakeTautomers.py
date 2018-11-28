import __future__

import random

import gypsum.parallelizer as parallelizer
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol
import gypsum.MolObjectHandling as MOH

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

try:
    from molvs import tautomer
except:
    Utils.log("You need to install molvs and its dependencies.")
    raise ImportError("You need to install molv and its dependencies.")

def make_tauts(contnrs, max_variants_per_compound, thoroughness, num_processors, multithread_mode, Parallelizer_obj):
    """
    Generates tautomers of the molecules. Note that some of the generated
    tautomers are not realistic. If you find a certain improbable
    substructure keeps popping up, add it to the list in the
    `prohibited_substructures` definition found with MyMol.py in the function remove_bizarre_substruc().
    """

    if max_variants_per_compound == 0:
        return

    Utils.log("Generating tautomers for all molecules...")

    params = []
    for contnr in contnrs:
        for mol_index, mol in enumerate(contnr.mols):
            params.append([contnr, mol_index, max_variants_per_compound])

    tmp = Parallelizer_obj.run(params, parallel_makeTaut, num_processors, multithread_mode)

    # Flatten the resulting list of lists
    #none_data = parallelizer.strip_none(tmp)
    none_data = tmp
    taut_data = parallelizer.flatten_list(none_data)

    # Remove bad tauts
    taut_data = tauts_no_break_arom_rngs(contnrs, taut_data, num_processors, multithread_mode, Parallelizer_obj)
    taut_data = tauts_no_elim_chiral(contnrs, taut_data, num_processors, multithread_mode, Parallelizer_obj)
    taut_data = tauts_no_change_hs_to_cs_unless_alpha_to_carbnyl(
        contnrs, taut_data, num_processors, multithread_mode, Parallelizer_obj
    )

    ChemUtils.bst_for_each_contnr_no_opt(contnrs, taut_data,
                                max_variants_per_compound, thoroughness)


def parallel_makeTaut(contnr, mol_index, max_variants_per_compound):
    #id = random.random()
    #print "Start", id

    mol = contnr.mols[mol_index]

    # Create a temporary mol, since that's what tauts works
    # with.
    # TODO: There should be a copy function
    m = MyMol.MyMol(mol.smiles()).rdkit_mol

    if m is None:
        Utils.log(
            "\tCould not generate tautomers for " + contnr.orig_smi +
            ". I'm deleting it."
        )
        return

    # Molecules should be kekulized already, but as sanity check
    # let's do that again. Because taut requires kekulized
    # input.
    Chem.Kekulize(m)
    m = MOH.check_sanitization(m)   
    if m is None:   
        return None   

    # Limit to max_variants_per_compound tauts. Note that another
    # batch could add more, so you'll need to once again trim to this
    # number later. But this could at least help prevent the
    # combinatorial explosion at this stage.
    enum = tautomer.TautomerEnumerator(
        max_tautomers=max_variants_per_compound
    )
    tauts_rdkit_mols = enum.enumerate(m)

    # Make all those tauts into MyMol objects
    tauts_mols = [MyMol.MyMol(m) for m in tauts_rdkit_mols]

    # Keep only those that have reasonable substructures
    tauts_mols = [t for t in tauts_mols if t.remove_bizarre_substruc() == False]

    if len(tauts_mols) > 1:
        Utils.log("\t" + mol.smiles(True) + " has tauts.")

    results = []

    for tm in tauts_mols:
        tm.inherit_contnr_props(contnr)
        tm.genealogy = mol.genealogy[:]
        tm.name = mol.name

        if tm.smiles() != mol.smiles():
            tm.genealogy.append(tm.smiles(True) + " (tautomer)")

        results.append(tm)

    return results
    #print "End", id


def tauts_no_break_arom_rngs(contnrs, taut_data, num_processors, multithread_mode, Parallelizer_obj):
    """
    For a given molecule, the number of atomatic rings should never change
    regardless of tautization, ionization, etc. Any taut that breaks
    aromaticity is unlikely to be worth pursuing. So remove it.

    :param [MyMol.MyMol] taut_data: A list of MyMol.MyMol objects.

    :returns: A list of MyMol.MyMol objects, with certain bad ones removed.
    :rtype: :class:`str` ???
    """

    # You need to group the taut_data by contnr
    params = []
    for taut_mol in taut_data:
        params.append([taut_mol, contnrs[taut_mol.contnr_idx]])

    tmp = Parallelizer_obj.run(params, parallel_CheckNonaroRings, 
                            num_processors, multithread_mode)

    # Stripping out None values
    results = parallelizer.strip_none(tmp)

    return results


def tauts_no_elim_chiral(contnrs, taut_data, num_processors, multithread_mode, Parallelizer_obj):
    """
    Unfortunately, molvs sees removing chiral specifications as being a
    distinct taut. I imagine there are cases where tautization could
    remove a chiral center, but I think these cases are rare. To compensate
    for the error in other folk's code, let's just require that the number of
    chiral centers remain unchanged with isomerization.

    :param [MyMol.MyMol] taut_data: A list of MyMol.MyMol objects.

    :returns: A list of MyMol.MyMol objects, with certain bad ones removed.
    :rtype: :class:`str` ???
    """

    # You need to group the taut_data by contnr
    params = []
    print("")
    print("")
    print("")
    print("contnrs: 165: ", contnrs)
    top_dir = "/ihome/jdurrant/jspiegel/gypsum/"
    import os
    import pickle
    if os.path.exists(top_dir) == False:
        top_dir = "/home/jacob/Documents/gypsum/"
        if os.path.exists(top_dir) == False:
            raise Exception("where is this?")
    pickle_file = top_dir + "picklefile"
    # if os.path.exists(pickle_file) == True:
    #     with open(pickle_file, 'rb') as handle:
    #             # should keep as list
    #             old_data = pickle.load(handle)
    with open(pickle_file, 'wb') as handle:
        pickle.dump(contnrs, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    with open(pickle_file+"taut_data", 'wb') as handle:
        pickle.dump(taut_data, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
    print("")
    print("")
    print("taut_data: ", taut_data)
    for taut_mol in taut_data:
        taut_mol_idx = int(taut_mol.contnr_idx)
        params.append([taut_mol, contnrs[taut_mol_idx]])
    print("")
    print("params: 165: ", params)
    print("")
    print("")
    print("")

    tmp = Parallelizer_obj.run(params, parallel_CheckChiralCenters, 
                            num_processors, multithread_mode)

    # Stripping out None values
    results = [x for x in tmp if x != None]

    return results

def tauts_no_change_hs_to_cs_unless_alpha_to_carbnyl(contnrs, taut_data, num_processors, multithread_mode, Parallelizer_obj):
    """
    Generally speaking, only carbons that are alpha to a carbonyl are
    sufficiently acidic to participate in taut formation. The
    taut-generating code you use makes these inappropriate tauts.

    :param [MyMol.MyMol] taut_data: A list of MyMol.MyMol objects.

    :returns: A list of MyMol.MyMol objects, with certain bad ones removed.
    :rtype: :class:`str` ???
    """

    # You need to group the taut_data by contnr
    params = []
    for taut_mol in taut_data:
        params.append([taut_mol, contnrs[taut_mol.contnr_idx]])

    tmp = Parallelizer_obj.run(params, parallel_CheckCarbonHydrogens, 
                            num_processors, multithread_mode)

    # Stripping out None values
    results = [x for x in tmp if x != None]

    return results



def parallel_CheckNonaroRings(taut, contnr):
    """
    A parallelizable helper function that checks that tautomers do not break any
    nonaromatic rings present in the original object.

    :param tautomer taut: A given tautomer to check.
    :param ?contnr? contnr: The original molecule object.

    :results: Returns either the tautomer or a None object.
    """

    # How many nonaromatic rings in the original smiles?
    num_nonaro_rngs_orig = contnr.num_nonaro_rngs

    # Check if it breaks aromaticity.
    m_num_nonaro_rngs = len(taut.m_num_nonaro_rngs())
    if m_num_nonaro_rngs == num_nonaro_rngs_orig:
        # Same number of nonaromatic rings as original molecule
        # Save the good ones.
        return taut
    else:
        Utils.log(
            "\t" + taut.smiles(True) + ", a tautomer generated " +
            "from " + contnr.orig_smi + " (" + taut.name +
            "), broke an aromatic ring, so I'm discarding it."
        )


def parallel_CheckChiralCenters(taut, contnr):
    """
    A parallelizable helper function that checks that tautomers do not break
    any chiral centers in the original molecule.

    :param tautomer taut: A given tautomer to check.
    :param ?contnr? contnr: The original molecule object.

    :results: Returns either the tautomer or a None object.
    """
    # How many chiral centers in the original smiles?
    num_specif_chiral_cntrs_orig = contnr.num_specif_chiral_cntrs

    # Make a new list containing only the ones that don't break chiral
    # centers.
    m_num_specif_chiral_cntrs = len(taut.chiral_cntrs_only_asignd())
    if m_num_specif_chiral_cntrs == num_specif_chiral_cntrs_orig:
        # Same number of chiral centers as original molecule
        # Save those good ones.
        return taut
    else:
        Utils.log(
            "\t" + contnr.orig_smi + " ==> " + taut.smiles(True) +
            " (tautomer transformation on " + taut.name + ") " +
            "changed the molecules total number of specified " +
            "chiral centers from " +
            str(num_specif_chiral_cntrs_orig) + " to " +
            str(m_num_specif_chiral_cntrs) +
            ", so I'm deleting it."
        )


def parallel_CheckCarbonHydrogens(taut, contnr):
    """
    A parallelizable helper function that checks that tautomers do not change
    the hydrogens on inappropriate carbons.

    :param tautomer taut: A given tautomer to check.
    :param ?contnr? contnr: The original molecule object.

    :results: Returns either the tautomer or a None object.
    """
    # What's the carbon-hydrogen fingerprint of the original smiles?
    orig_carbon_hydrogen_count = contnr.carbon_hydrogen_count

    # How about this taut?
    this_carbon_hydrogen_count = taut.carb_hyd_cnt()

    # Only keep if they are the same.
    if orig_carbon_hydrogen_count == this_carbon_hydrogen_count:
        return taut
    else:
        Utils.log(
            "\t" + contnr.orig_smi + " ==> " + taut.smiles(True) +
            " (taut transformation on " + taut.name + ") " +
            "changed the number of hydrogen atoms bound to a " +
            "carbon, so I'm deleting it."
        )
