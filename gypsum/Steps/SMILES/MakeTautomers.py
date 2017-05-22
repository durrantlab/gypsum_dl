import sys
import random

import gypsum.Multiprocess as mp
import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    sys.exit(0)

try:
    from molvs import tautomer
except:
    Utils.log("You need to install molvs and its dependencies.")
    sys.exit(0)

def make_tauts(self):
    """
    Generates tautomers of the molecules. Note that some of the generated
    tautomers are not realistic. If you find a certain improbable
    substructure keeps popping up, add it to the list in the
    `remove_crazy_substrings` definition.
    """

    if self.params["max_variants_per_compound"] == 0:
        return

    Utils.log("Generating tautomers for all molecules...")

    params = []
    for contnr in self.contnrs:
        for mol_index, mol in enumerate(contnr.mols):
            params.append((contnr, mol_index, self.params))

    tmp = mp.MultiThreading(params, self.params["num_processors"], parallel_makeTaut)

    # Flatten the resulting list of lists
    #none_data = mp.strip_none(tmp)
    none_data = tmp
    print tmp
    taut_data = mp.flatten_list(none_data)
    print taut_data

    # Remove bad tauts
    taut_data = tauts_no_break_arom_rngs(self, taut_data)
    taut_data = tauts_no_elim_chiral(self, taut_data)
    taut_data = tauts_no_change_hs_to_cs_unless_alpha_to_carbnyl(
        self, taut_data
    )
    print taut_data

    ChemUtils.bst_for_each_contnr_no_opt(self, taut_data)


def parallel_makeTaut(contnr, mol_index, params):
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

    # Limit to max_variants_per_compound tauts. Note that another
    # batch could add more, so you'll need to once again trim to this
    # number later. But this could at least help prevent the
    # combinatorial explosion at this stage.
    enum = tautomer.TautomerEnumerator(
        max_tautomers=params["max_variants_per_compound"]
    )
    tauts_rdkit_mols = enum.enumerate(m)

    # Make all those tauts into MyMol objects
    tauts_mols = [MyMol.MyMol(m) for m in tauts_rdkit_mols]

    # Keep only those that have reasonable substructures
    tauts_mols = [t for t in tauts_mols if t.crzy_substruc() == False]

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


def tauts_no_break_arom_rngs(self, taut_data):
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
        params.append((taut_mol, self.contnrs[taut_mol.contnr_idx]))

    tmp = mp.MultiThreading(params, self.params["num_processors"],
                            parallel_CheckNonaroRings)

    # Stripping out None values
    results = mp.strip_none(tmp)

    return results


def tauts_no_elim_chiral(self, taut_data):
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
    for taut_mol in taut_data:
        params.append((taut_mol, self.contnrs[taut_mol.contnr_idx]))

    tmp = mp.MultiThreading(params, self.params["num_processors"],
                            parallel_CheckChiralCenters)

    # Stripping out None values
    results = [x for x in tmp if x != None]

    return results

def tauts_no_change_hs_to_cs_unless_alpha_to_carbnyl(self, taut_data):
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
        params.append((taut_mol, self.contnrs[taut_mol.contnr_idx]))

    tmp = mp.MultiThreading(params, self.params["num_processors"],
                            parallel_CheckCarbonHydrogens)

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
