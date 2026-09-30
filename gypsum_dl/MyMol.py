"""
This module contains classes and functions for processing individual molecules
(variants). All variants of the same input molecule are grouped together in
the same MolContainer.MolContainer object. Each MyMol.MyMol is also associated
with conformers described here (3D coordinate sets).

So just to clarify: MolContainer.MolContainer > MyMol.MyMol >
MyMol.MyConformer
"""

import contextlib
import copy
import operator
import random
import sys
from typing import TYPE_CHECKING, TypedDict

from molvs import standardize_smiles as ssmiles

# Disable the unnecessary RDKit warnings
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondStereo

import gypsum_dl.MolObjectHandling as MOH
from gypsum_dl import utils

RDLogger.DisableLog("rdApp.*")

if TYPE_CHECKING:
    # Annotations only; MolContainer imports this module.
    from gypsum_dl.MolContainer import MolContainer


class ContnrProps(TypedDict):
    """The container-level fields every variant of a compound shares.

    These fields describe the input compound, not the individual variant, so a
    step that needs them does not need the container itself. Naming them in one
    place is what lets a job carry them: the tautomer step used to ship a whole
    MolContainer (and so every variant in it) as an argument to each per-variant
    job, which made its pickling cost grow with the square of the variant count.
    Having a single definition also keeps the two ways a variant acquires these
    fields, MolContainer.add_smiles and inherit_contnr_props, from drifting
    apart on which fields they set.
    """

    # A MyMol that no container has claimed yet carries "" here, so the union
    # describes what the field actually holds rather than what it ought to.
    contnr_idx: int | str
    name: str
    orig_smi: str
    orig_smi_deslt: str
    orig_smi_canonical: str | None


class _Unset:
    """Marker for a memoized field that has not been computed yet.

    Several of the accessors below cache values that are themselves falsy or
    None: remove_bizarre_substruc caches False, get_idxs_of_nonaro_rng_atms
    caches an empty list, and smiles() caches None to record a failed
    canonicalization. An empty string used to serve as "not computed yet",
    which left every cached field with three states and made each read site
    re-derive which falsy value meant what. A marker that is none of those
    values reduces the question to `is UNSET`.
    """

    def __bool__(self) -> bool:
        """Refuse to act as a truth value.

        `if not self.frgs:` is the wrong question to ask of one of these
        fields, since a legitimately cached False or empty list answers it the
        same way an unfilled cache does. Raising turns that mistake into an
        immediate error rather than a cache that is silently recomputed or
        silently returned as a result.

        Raises:
            TypeError: Always.
        """

        raise TypeError(
            "An unfilled Gypsum-DL cache has no truth value; compare it with "
            "`is UNSET` instead."
        )

    def __repr__(self) -> str:
        """Name the marker in tracebacks and log output.

        Returns:
            The name this module exports the marker under.
        """

        return "UNSET"

    def __reduce__(self) -> str:
        """Stay a singleton across pickling and copying.

        MyMol objects are pickled to and from worker processes at every stage
        of the pipeline, and deep-copied besides. Reconstructing this marker as
        a separate instance would make `is UNSET` report every unfilled cache
        as filled on the far side of the parallelizer. Returning the global
        name has pickle and copy resolve back to this instance.

        Returns:
            The module-level name bound to this instance.
        """

        return "UNSET"


UNSET = _Unset()


def _atom_order_preserved(before: Chem.Mol, after: Chem.Mol) -> bool:
    """Report whether every atom of `before` keeps its index in `after`.

    MyMol caches atom indices (ring members, chiral centers), which remain
    meaningful across a molecule swap only if the old atoms are a prefix of
    the new ones.

    Args:
        before: The molecule the caches were computed on.
        after: The molecule replacing it.

    Returns:
        True if `after` starts with the atoms of `before`, in order, with the
        same elements and with every bond of `before` joining the same indices.
    """

    if after.GetNumAtoms() < before.GetNumAtoms():
        return False
    same_elements = all(
        after.GetAtomWithIdx(i).GetAtomicNum() == atom.GetAtomicNum()
        for i, atom in enumerate(before.GetAtoms())
    )
    # Elements alone would miss a permutation among atoms of the same kind.
    return same_elements and all(
        after.GetBondBetweenAtoms(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) is not None
        for b in before.GetBonds()
    )


class MyMol:
    """
    A class that wraps around a rdkit.Mol object. Includes additional data and
    functions.
    """

    def __init__(self, starter, name=""):
        """Initialize the MyMol object.

        :param starter: The object (smiles or rdkit.Mol) on which to build this
           class.
        :type starter: str or rdkit.Mol
        :param name: An optional string, the name of this molecule. Defaults to "".
        :param name: str, optional
        """

        if isinstance(starter, str):
            # It's a SMILES string.
            self.rdkit_mol = ""
            self.can_smi = UNSET
            smiles = starter
        else:
            # So it's an rdkit mol object. No need to regenerate it, but do
            # take a copy: this object mutates its molecule in place
            # (sanitization, reprotonation, and above all
            # load_conformers_into_rdkit_mol, which clears and rewrites the
            # conformer set), so holding the caller's object would let those
            # edits appear underneath whoever passed it in. The isinstance
            # guard keeps the graceful path for a starter that is neither a
            # string nor a molecule: MakeTautomers builds a MyMol from
            # smiles(), which is None once canonicalization has failed.
            self.rdkit_mol = (
                Chem.Mol(starter) if isinstance(starter, Chem.Mol) else starter
            )

            # Get the smiles too from the rdkit mol object.
            try:
                smiles = Chem.MolToSmiles(
                    self.rdkit_mol, isomericSmiles=True, canonical=True
                )

                # In this case you know it's cannonical.
                self.can_smi = smiles
            except Exception:
                # Sometimes this conversion just can't happen. Happened once
                # with this beast, for example:
                # CC(=O)NC1=CC(=C=[N+]([O-])O)C=C1O
                # None is the same failure marker smiles() uses, so callers
                # have one non-string case to reason about rather than two.
                self.can_smi = None
                smiles = ""
                id_to_print = name if name != "" else str(starter)
                utils.log(
                    "\tERROR: Could not generate one of the structures "
                    + "for ("
                    + id_to_print
                    + ")."
                )

        self.can_smi_noh = UNSET
        self.orig_smi = smiles

        # Default assumption is that they are the same.
        self.orig_smi_deslt = smiles

        # A container-level field, so it is None until a container claims this
        # molecule (MolContainer.add_smiles and inherit_contnr_props both fill
        # it). Declared here because it was previously set only by add_smiles,
        # which left which fields exist on a MyMol depending on which code path
        # built it.
        self.orig_smi_canonical: str | None = None
        self.name = name
        self.conformers = []
        self.nonaro_ring_atom_idx = UNSET
        self.chiral_cntrs_only_assigned = UNSET
        self.chiral_cntrs_include_unasignd = UNSET
        self.bizarre_substruct = UNSET
        self.enrgy = {}  # different energies for different conformers.
        self.minimized_enrgy = {}
        self.contnr_idx = ""
        self.frgs = UNSET
        self.stdrd_smiles = UNSET
        self.mol_props = {}
        self.idxs_low_energy_confs_no_opt = {}
        self.idxs_of_confs_to_min = set([])
        self.genealogy = []  # Keep track of how the molecule came to be.

        # Makes the molecule if a smiles was provided. Sanitizes the molecule
        # regardless.
        sanitized = self.make_mol_frm_smiles_sanitze()

        if isinstance(self.can_smi, str):
            # check_sanitization can hand back a modified copy of the molecule
            # (the four-bond nitrogen fix), and sanitizing in place can change
            # how a molecule is written, so any SMILES canonicalized above may
            # describe a molecule this object no longer holds. Drop it so
            # smiles() recomputes from self.rdkit_mol. If sanitization failed
            # outright there is nothing to recompute from (rdkit_mol is None),
            # so record the failure the way smiles() does rather than keeping a
            # string that describes the structure we just discarded: callers
            # that guard on a non-string SMILES (uniq_mols_in_list,
            # contains_canonical_smiles, remove_highly_charged_molecules) would
            # otherwise treat this object as a usable molecule.
            self.can_smi = UNSET if sanitized is not None else None

    def standardize_smiles(self):
        """Standardize the smiles string if you can."""

        if self.stdrd_smiles is not UNSET:
            return self.stdrd_smiles

        try:
            self.stdrd_smiles = ssmiles(self.smiles())
        except Exception:
            # orig_smi is always a string, whereas smiles(True) reports a
            # failed deprotonation as None, so building the message from it
            # replaced the original failure with a TypeError raised inside
            # this handler.
            utils.log(
                f"\tCould not standardize {self.orig_smi} ({self.name}). Skipping."
            )
            self.stdrd_smiles = self.smiles()

        return self.stdrd_smiles

    def __hash__(self):
        """Allows you to compare MyMol.MyMol objects.

        :return: The hashed canonical smiles.
        :rtype: str
        """

        can_smi = self.smiles()

        if not isinstance(can_smi, str):
            # smiles() reports failure as None. Hashing that would give every
            # molecule whose canonical SMILES could not be computed the same
            # hash, so the deduplication passes (dict.fromkeys here,
            # uniq_mols_in_list elsewhere) would treat unrelated molecules as
            # copies of each other and keep only one. Fall back to identity so
            # such a molecule is equal only to itself.
            return object.__hash__(self)

        # So it hashes based on the cannonical smiles.
        return hash(can_smi)

    def sort_key(self) -> tuple[int, str]:
        """Give a stable key for ordering molecules.

        The ordering operators used to compare hashes, but str hashing is
        salted per interpreter, so any sort that fell through to comparing
        molecules (sorting (energy, MyMol) pairs whose energies tie, for
        instance) put them in a different order on every run. The canonical
        SMILES is stable across runs. A molecule that cannot be canonicalized
        sorts after every molecule that can, and ties with the others that
        cannot, which leaves those in their original order under a stable
        sort.

        Returns:
            A tuple of (rank, canonical smiles), comparable against the same
                from any other molecule.
        """

        can_smi = self.smiles()

        if not isinstance(can_smi, str):
            return (1, "")

        return (0, can_smi)

    def __eq__(self, other):
        """Allows you to compare MyMol.MyMol objects.

        :param other: The other molecule.
        :type other: MyMol.MyMol
        :return: Whether the other molecule is the same as this one.
        :rtype: bool
        """

        if not isinstance(other, MyMol):
            # Anything that is not a molecule is not this molecule, and asking
            # it for a hash is misleading besides: hash(True) is 1, so a
            # molecule could otherwise compare equal to a bool sentinel.
            return False

        can_smi = self.smiles()
        other_can_smi = other.smiles()

        if not isinstance(can_smi, str) or not isinstance(other_can_smi, str):
            # smiles() reports failure as None. Two molecules that both failed
            # to canonicalize have not been shown to be the same molecule, so
            # such a molecule is equal only to itself. This matches the
            # identity fallback in __hash__.
            return self is other

        # Compare the canonical smiles themselves rather than their hashes:
        # two distinct molecules whose SMILES happen to hash to the same value
        # are not the same molecule, and the deduplication passes treat
        # equality as "this is a copy, drop it."
        return can_smi == other_can_smi

    def __ne__(self, other):
        """Allows you to compare MyMol.MyMol objects.

        :param other: The other molecule.
        :type other: MyMol.MyMol
        :return: Whether the other molecule is different from this one.
        :rtype: bool
        """

        return not self.__eq__(other)

    def __lt__(self, other):
        """Is this MyMol less than another one? Gypsum-DL often sorts
        molecules by sorting tuples of the form (energy, MyMol). On rare
        occasions, the energies are identical, and the sorting algorithm
        attempts to compare MyMol directly.

        :param other: The other molecule.
        :type other: MyMol.MyMol
        :return: True or False, if less than or not.
        :rtype: boolean
        """

        if not isinstance(other, MyMol):
            return NotImplemented

        return self.sort_key() < other.sort_key()

    def __le__(self, other):
        """Is this MyMol less than or equal to another one? Gypsum-DL often
        sorts molecules by sorting tuples of the form (energy, MyMol). On rare
        occasions, the energies are identical, and the sorting algorithm
        attempts to compare MyMol directly.

        :param other: The other molecule.
        :type other: MyMol.MyMol
        :return: True or False, if less than or equal to, or not.
        :rtype: boolean
        """

        if not isinstance(other, MyMol):
            return NotImplemented

        # Deferring to the strict operator and equality keeps the two in step:
        # molecules that merely tie in the sort order (two that could not be
        # canonicalized) are not equal, so neither is <= the other.
        return self.__lt__(other) or self.__eq__(other)

    def __gt__(self, other):
        """Is this MyMol greater than another one? Gypsum-DL often sorts
        molecules by sorting tuples of the form (energy, MyMol). On rare
        occasions, the energies are identical, and the sorting algorithm
        attempts to compare MyMol directly.

        :param other: The other molecule.
        :type other: MyMol.MyMol
        :return: True or False, if greater than or not.
        :rtype: boolean
        """

        if not isinstance(other, MyMol):
            return NotImplemented

        return self.sort_key() > other.sort_key()

    def __ge__(self, other):
        """Is this MyMol greater than or equal to another one? Gypsum-DL often
        sorts molecules by sorting tuples of the form (energy, MyMol). On rare
        occasions, the energies are identical, and the sorting algorithm
        attempts to compare MyMol directly.

        :param other: The other molecule.
        :type other: MyMol.MyMol
        :return: True or False, if greater than or equal to, or not.
        :rtype: boolean
        """

        if not isinstance(other, MyMol):
            return NotImplemented

        return self.__gt__(other) or self.__eq__(other)

    def make_mol_frm_smiles_sanitze(self):
        """Construct a rdkit.mol for this object, in case you only received
        the smiles. Also, sanitize the molecule regardless.

        :return: Returns the rdkit.mol object, though it's also stored in
           self.rdkit_mol.
        :rtype: rdkit.mol object.
        """

        # If given a SMILES string.
        if self.rdkit_mol == "":
            try:
                # sanitize = False makes it respect double-bond stereochemistry
                m = Chem.MolFromSmiles(self.orig_smi_deslt, sanitize=False)
            except Exception:
                m = None
        else:  # If given a RDKit Mol Obj
            m = self.rdkit_mol

        if m is not None:
            # Sanitize and hopefully correct errors in the smiles string such
            # as incorrect nitrogen charges.
            m = MOH.check_sanitization(m)
        self.rdkit_mol = m
        return m

    def make_first_3d_conf_no_min(self):
        """Makes the associated rdkit.mol object 3D by adding the first
        conformer. This also adds hydrogen atoms to the associated rdkit.mol
        object. Note that it does not perform a minimization, so it is not
        too expensive."""

        # Set the first 3D conformer
        if len(self.conformers) > 0:
            # It's already been done.
            return

        # Add hydrogens. This adds explicit hydrogens, while respecting
        # Dimorphite-DL protonation states. Reprotonation can fail and return
        # None; assigning that to rdkit_mol would crash MyConformer downstream
        # (deepcopy(None).RemoveAllConformers()), so bail out quietly instead.
        reprotanated = MOH.try_reprotanation(self.rdkit_mol)
        if reprotanated is None:
            return

        # Fill the canonical-SMILES cache before swapping in the reprotonated
        # molecule. AddHs puts the hydrogens into the graph, so canonicalizing
        # afterwards writes them out explicitly ("[H]O[H]" rather than "O"),
        # and whether that happened depended on whether anything had called
        # smiles() earlier: add_conformers below reaches MyConformer, which
        # reads smiles() and so fills an empty cache from the hydrogen-added
        # molecule. The cached value is this molecule's identity key
        # (deduplication, the SMILES recorded in the output), and adding
        # explicit hydrogens does not change its identity, so pin it to the
        # form the rest of the pipeline already uses.
        self.smiles()

        # The ring and chiral-center caches hold atom indices. AddHs appends
        # the hydrogens after the existing atoms, so those indices stay valid
        # and keeping the caches leaves results unchanged; anything that
        # renumbers the heavy atoms would silently invalidate them. The
        # bizarre-substructure verdict names no atoms, so like the SMILES
        # above it stays pinned to the molecule as it was before the swap.
        if not _atom_order_preserved(self.rdkit_mol, reprotanated):
            self.nonaro_ring_atom_idx = UNSET
            self.chiral_cntrs_only_assigned = UNSET
            self.chiral_cntrs_include_unasignd = UNSET

        self.rdkit_mol = reprotanated

        # Add a single conformer. RMSD cutoff very small so all conformers
        # will be accepted. And not minimizing (False).
        self.add_conformers(1, 1e60, False)

    def smiles(self, noh=False):
        """Get the desalted, canonical smiles string associated with this
           object. (Not the input smiles!)

        :param noh: Whether or not hydrogen atoms should be included in the
           canonical smiles string., defaults to False
        :param noh: bool, optional
        :return: The canonical smiles string, or None if it cannot be
           determined.
        :rtype: str or None
        """

        # See if it's already been calculated. They want the hydrogen atoms.
        if noh == False:
            if self.can_smi is not UNSET:
                # Return previously determined canonical SMILES.
                return self.can_smi

            # Need to determine canonical SMILES.
            try:
                can_smi = Chem.MolToSmiles(
                    self.rdkit_mol, isomericSmiles=True, canonical=True
                )
            except Exception:
                utils.log(
                    f"Warning: Couldn't put {self.orig_smi} ({self.name}) in canonical form. Got this error: {str(sys.exc_info()[0])}. This molecule will be discarded."
                )
                self.can_smi = None
                return None

            self.can_smi = can_smi
            return can_smi
        else:
            # They don't want the hydrogen atoms.
            if self.can_smi_noh is not UNSET:
                # Return previously determined string.
                return self.can_smi_noh

            # So remove hydrogens. Note that this assumes you will have called
            # this function previously with noh = False
            amol = copy.copy(self.rdkit_mol)
            amol = MOH.try_deprotanation(amol)
            if amol is None:
                # rdkit_mol can be None (a molecule that failed sanitization or
                # 3D optimization), and deprotonation can fail on its own.
                # MolToSmiles would then raise, and unlike the noh == False
                # branch, several callers of this accessor run in the main
                # process after every expensive step, so the exception aborts
                # the whole run before any output is written. Report failure the
                # same way the other branch does instead.
                utils.log(
                    f"Warning: Couldn't put {self.orig_smi} ({self.name}) in canonical form without hydrogens. This molecule will be discarded."
                )
                self.can_smi_noh = None
                return None
            self.can_smi_noh = Chem.MolToSmiles(
                amol, isomericSmiles=True, canonical=True
            )
            return self.can_smi_noh

    def get_idxs_of_nonaro_rng_atms(self):
        """Identifies which rings in a given molecule are nonaromatic, if any.

        :return: A [[int, int, int]]. A list of lists, where each inner list is
           a list of the atom indecies of the members of a non-aromatic ring.
           Also saved to self.nonaro_ring_atom_idx.
        :rtype: list
        """

        if self.nonaro_ring_atom_idx is not UNSET:
            # Already determined...
            return self.nonaro_ring_atom_idx

        # There are no rings if the molecule is None.
        if self.rdkit_mol is None:
            return []

        # Get the number of symmetric smallest set of rings
        ssr = Chem.GetSymmSSSR(self.rdkit_mol)

        # Get the rings
        ring_indecies = [list(ssr[i]) for i in range(len(ssr))]

        # Are the atoms in any of those rings nonaromatic?
        nonaro_rngs = []
        for rng_indx_set in ring_indecies:
            for atm_idx in rng_indx_set:
                if self.rdkit_mol.GetAtomWithIdx(atm_idx).GetIsAromatic() == False:
                    # One of the ring atoms is not aromatic! Let's keep it.
                    nonaro_rngs.append(rng_indx_set)
                    break
        self.nonaro_ring_atom_idx = nonaro_rngs
        return nonaro_rngs

    def chiral_cntrs_w_unasignd(self):
        """Get every chiral center, whether or not it has been assigned.

        Despite the name, this is a superset of chiral_cntrs_only_asignd()
        rather than its complement: unassigned centers are marked '?' and
        assigned ones carry their 'R' or 'S'. So the length of this list is the
        total number of chiral centers, and adding it to the length of
        chiral_cntrs_only_asignd() counts every assigned center twice. Callers
        wanting only the unassigned centers filter on '?' (see
        EnumerateChiralMols.parallel_get_chiral).

        :return: The chiral centers. Also saved to
           self.chiral_cntrs_include_unasignd. Looks like [(10, '?')]
        :rtype: list
        """

        # No chiral centers if the molecule is None.
        if self.rdkit_mol is None:
            return []

        if self.chiral_cntrs_include_unasignd is not UNSET:
            # Already been determined...
            return self.chiral_cntrs_include_unasignd

        # Get the chiral centers that are not defined. The perception mode
        # decides which atoms count, and the tautomer filter compares these
        # counts across molecules, so it must not float with RDKit's default.
        with MOH.legacy_stereo_perception():
            ccs = Chem.FindMolChiralCenters(self.rdkit_mol, includeUnassigned=True)
        self.chiral_cntrs_include_unasignd = ccs
        return ccs

    def chiral_cntrs_only_asignd(self):
        """Get the chiral centers that have been assigned.

        :return: The chiral centers. Also saved to self.chiral_cntrs_only_assigned.
        :rtype: list
        """

        if self.chiral_cntrs_only_assigned is not UNSET:
            return self.chiral_cntrs_only_assigned

        if self.rdkit_mol is None:
            return []

        with MOH.legacy_stereo_perception():
            ccs = Chem.FindMolChiralCenters(self.rdkit_mol, includeUnassigned=False)
        self.chiral_cntrs_only_assigned = ccs
        return ccs

    def get_double_bonds_without_stereochemistry(self):
        """Get the double bonds that don't have specified stereochemistry.

        :return: The unasignd double bonds (indexes). Looks like this:
           [2, 4, 7]
        :rtype: list
        """

        if self.rdkit_mol is None:
            return []

        return [
            b.GetIdx()
            for b in self.rdkit_mol.GetBonds()
            if b.GetBondTypeAsDouble() == 2 and b.GetStereo() is BondStereo.STEREONONE
        ]

    def remove_bizarre_substruc(self):
        """Removes molecules with improbable substuctures, likely generated
           from the tautomerization process. Used to find artifacts.

        :return: Boolean, whether or not there are impossible substructures.
           Also saves to self.bizarre_substruct.
        :rtype: bool
        """

        if self.bizarre_substruct is not UNSET:
            # Already been determined.
            return self.bizarre_substruct

        if self.rdkit_mol is None:
            # It is bizarre to have a molecule with no atoms in it.
            return True

        # These are substrutures that can't be easily corrected using
        # fix_common_errors() below.
        # , "[C+]", "[C-]", "[c+]", "[c-]", "[n-]", "[N-]"] # ,
        # "[*@@H]1(~[*][*]~2)~[*]~[*]~[*@@H]2~[*]~[*]~1",
        # "[*@@H]1~2~*~*~[*@@H](~*~*2)~*1",
        # "[*@@H]1~2~*~*~*~[*@@H](~*~*2)~*1",
        # "[*@@H]1~2~*~*~*~*~[*@@H](~*~*2)~*1",
        # "[*@@H]1~2~*~[*@@H](~*~*2)~*1", "[*@@H]~1~2~*~*~*~[*@H]1O2",
        # "[*@@H]~1~2~*~*~*~*~[*@H]1O2"]

        # Note that C(O)=N, C and N mean they are aliphatic. Does not match
        # c(O)n, when aromatic. So this form is acceptable if in aromatic
        # structure.
        prohibited_substructures = ["O(=*)-*"]  # , "C(O)=N"]

        # Enol forms with terminal alkenes are unlikely.
        prohibited_substructures.append("C(=[CH2])[OH]")

        # Enol forms with terminal alkenes are unlikely.
        prohibited_substructures.append("C(=[CH2])[O-]")
        # A geminal vinyl diol is not a tautomer of a carboxylate group.
        prohibited_substructures.append("C=C([OH])[OH]")

        # A geminal vinyl diol is not a tautomer of a carboxylate group.
        prohibited_substructures.append("C=C([O-])[OH]")

        # A geminal vinyl diol is not a tautomer of a carboxylate group.
        prohibited_substructures.append("C=C([O-])[O-]")
        prohibited_substructures.append("[C-]")  # No carbanions.
        prohibited_substructures.append("[c-]")  # No carbanions.

        # The patterns above are SMARTS, so they are matched against the
        # molecule and not against its SMILES strings. A substring pass used
        # to run first, on the theory that it was a cheap approximation, but
        # most of these patterns cannot appear verbatim in SMILES output
        # (RDKit writes [OH] as O, and [C-] is not a substring of [CH2-]), and
        # the two that can are matched here anyway. Leaving it in place
        # suggested that a pattern added to the list would be enforced by
        # string comparison, which it would not be.
        for s in prohibited_substructures:
            pattrn = Chem.MolFromSmarts(s)
            if self.rdkit_mol.HasSubstructMatch(pattrn):
                # utils.log("\tRemoving a molecule because it has an odd
                # substructure: " + s)
                utils.log("\tDetected unusual substructure: " + s)
                self.bizarre_substruct = True
                return True

        # Now certin patterns that are more complex.
        # TODO in the future?

        self.bizarre_substruct = False
        return False

    def get_frags_of_orig_smi(self) -> list[Chem.Mol]:
        """Divide the current molecule into fragments.

        The single-fragment branch used to return [self], so the element type
        depended on the fragment count and the only caller (the desalter) was
        safe only because that branch returns before anything reads an rdkit
        method off the elements. Returning the wrapped molecule keeps the list
        homogeneous, and both the fragment count and GetNumHeavyAtoms work on
        it either way.

        Returns:
            The fragments as rdkit mols. Also caches them on self.frgs.
        """

        if self.frgs is not UNSET:
            # Already been determined...
            return self.frgs

        if "." not in self.orig_smi:
            # There are no fragments. Just return this molecule.
            self.frgs = [self.rdkit_mol]
            return self.frgs

        # Get the fragments.
        self.frgs = list(Chem.GetMolFrags(self.rdkit_mol, asMols=True))
        return self.frgs

    def contnr_props(self) -> ContnrProps:
        """Collect the container-level fields this molecule carries.

        Lets a step hand a job the few fields that describe the input compound
        instead of the MolContainer holding every variant of it. See
        ContnrProps.

        Returns:
            This molecule's copy of the container-level fields.
        """

        return {
            "contnr_idx": self.contnr_idx,
            "name": self.name,
            "orig_smi": self.orig_smi,
            "orig_smi_deslt": self.orig_smi_deslt,
            "orig_smi_canonical": self.orig_smi_canonical,
        }

    def inherit_contnr_props(self, other: "MolContainer | MyMol | ContnrProps") -> None:
        """Copies a few key properties from a different MyMol.MyMol object to
           this one.

        :param other: The container, molecule, or already-extracted mapping of
           container-level fields to copy from.
        :type other: MolContainer.MolContainer | MyMol.MyMol | ContnrProps
        """

        # These are properties that should be the same for every MyMol.MyMol
        # object in this MolContainer. A mapping is accepted so a step that
        # ships these fields to a worker in place of the container can still
        # apply them on the far side, without a second copy of the field list.
        props = other if isinstance(other, dict) else other.contnr_props()

        self.contnr_idx = props["contnr_idx"]
        self.orig_smi = props["orig_smi"]
        self.orig_smi_deslt = props["orig_smi_deslt"]  # initial assumption
        self.orig_smi_canonical = props["orig_smi_canonical"]
        self.name = props["name"]

    def set_rdkit_mol_prop(self, key: str, val: object) -> None:
        """Set a molecular property.

        :param key: The name of the molecular property.
        :type key: str
        :param val: The value of that property. A value of None is skipped
           rather than written, so a failed calculation (smiles(True) returning
           None, say) leaves the field absent instead of writing the literal
           string "None" into the output.
        :type val: object
        """

        if val is None:
            return

        val = str(val)
        with contextlib.suppress(Exception):
            self.rdkit_mol.SetProp(key, val)

    def set_all_rdkit_mol_props(self):
        """Set all the stored molecular properties. Copies ones from the
        MyMol.MyMol object to the MyMol.rdkit_mol object."""

        # self.set_rdkit_mol_prop("SOURCE_SMILES", self.orig_smi)
        for prop in list(self.mol_props.keys()):
            self.set_rdkit_mol_prop(prop, self.mol_props[prop])

        # SMILES, Genealogy, and _Name describe this variant, so they are
        # written after mol_props. An input SDF can carry a tag named SMILES,
        # and mol_props holds those input tags; writing SMILES first let the
        # input value replace the prepared variant's own SMILES.
        self.set_rdkit_mol_prop("SMILES", self.smiles(True))
        genealogy = "\n".join(self.genealogy)
        self.set_rdkit_mol_prop("Genealogy", genealogy)
        self.set_rdkit_mol_prop("_Name", self.name)

    def add_conformers(
        self,
        num: int,
        rmsd_cutoff: float = 0.1,
        minimize: bool = True,
        second_embed: bool = False,
    ) -> None:
        """Add conformers to this molecule.

        :param num: The total number of conformers to generate, including ones
           that have been generated previously.
        :type num: int
        :param rmsd_cutoff: Don't keep conformers that come within this rms
           distance of other conformers. Defaults to 0.1
        :param rmsd_cutoff: float, optional
        :param minimize: Whether or not to minimize the geometry of all these
           conformers. Defaults to True.
        :param minimize: bool, optional
        :param second_embed: Whether each new conformer may fall back to a
           last-resort embedding attempt if the default ones fail. Defaults to
           False.
        :type second_embed: bool, optional
        """

        # First, do you need to add new conformers? Some might have already
        # been added. Just add enough to meet the requested amount.
        num_new_confs = max(0, num - len(self.conformers))
        for _ in range(num_new_confs):
            if len(self.conformers) == 0:
                # For the first one, don't start from random coordinates.
                new_conf = MyConformer(self, None, second_embed)
            else:
                # For all subsequent ones, do start from random coordinates.
                new_conf = MyConformer(self, None, second_embed, True)

            if new_conf.mol is not False:
                self.conformers.append(new_conf)

        # Are the current ones minimized if necessary?
        if minimize == True:
            for conf in self.conformers:
                conf.minimize()  # Won't reminimize if it's already been done.

        # Automatically sort by the energy.
        self.conformers.sort(key=operator.attrgetter("energy"))

        # print(len(self.conformers))

        # Print the coordinates of the atoms of each conformer
        # for conf in self.conformers:
        #     print(conf.coords())
        #     print("")

        # Save all the conformers to separate PDB files
        # for i, conf in enumerate(self.conformers):
        #     conf.write_pdb_file("test_" + str(i) + ".pdb")
        #     print(conf.get_energy())

        # Remove ones that are very structurally similar.
        self.eliminate_structurally_similar_conformers(rmsd_cutoff)

    def eliminate_structurally_similar_conformers(self, rmsd_cutoff=0.1):
        """Eliminates conformers that are very geometrically similar.

        :param rmsd_cutoff: The RMSD cutoff to use. Defaults to 0.1
        :param rmsd_cutoff: float, optional
        """

        # Eliminate redundant ones.
        for i1 in range(len(self.conformers) - 1):
            if self.conformers[i1] is not None:
                for i2 in range(i1 + 1, len(self.conformers)):
                    if self.conformers[i2] is not None:
                        # Align them.
                        self.conformers[i2] = self.conformers[i1].align_to_me(
                            self.conformers[i2]
                        )

                        # Calculate the RMSD.
                        rmsd = self.conformers[i1].rmsd_to_me(self.conformers[i2])

                        # Replace the second one with None if it's too similar
                        # to the first.
                        if rmsd <= rmsd_cutoff:
                            self.conformers[i2] = None

        # Remove all the None entries.
        while None in self.conformers:
            self.conformers.remove(None)

        # Those that remains are only the distinct conformers.

    def load_conformers_into_rdkit_mol(self):
        """Load the conformers stored as MyConformers objects (in
        self.conformers) into the rdkit Mol object."""

        if self.rdkit_mol is None:
            return
        self.rdkit_mol.RemoveAllConformers()
        if not self.conformers:
            # Without this, SDWriter emits a coordinate block of zeros, so a
            # file advertised as 2D output would read as a degenerate 3D
            # structure with all atoms stacked at the origin.
            AllChem.Compute2DCoords(self.rdkit_mol)
            return
        for conformer in self.conformers:
            # Freshly embedded conformers all carry id 0, so without assignId
            # a multi-conformer molecule ends up with several conformers
            # sharing an id, and writers that resolve confId=-1 to the first
            # match would silently write the same coordinates repeatedly.
            self.rdkit_mol.AddConformer(conformer.conformer(), assignId=True)


class MyConformer:
    """A wrapper around a rdkit Conformer object. Allows me to associate extra
    values with conformers. These are 3D coordinate sets for a given
    MyMol.MyMol object (different molecule conformations).
    """

    def __init__(
        self, mol, conformer=None, second_embed=False, use_random_coordinates=False
    ):
        """Create a MyConformer objects.

        :param mol: The MyMol.MyMol associated with this conformer.
        :type mol: MyMol.MyMol
        :param conformer: An optional variable specifying the conformer to use.
           If not specified, it will create a new conformer. Defaults to None.
        :type conformer: rdkit.Conformer, optional
        :param second_embed: Whether to try to generate 3D coordinates using an
            older algorithm if the better (default) algorithm fails. This can add
            run time, but sometimes converts certain molecules that would
            otherwise fail. Defaults to False.
        :type second_embed: bool, optional
        :param use_random_coordinates: The first conformer should not start
           from random coordinates, but rather the eigenvalues-based
           coordinates rdkit defaults to. But Gypsum-DL generates subsequent
           conformers to try to consider alternate geometries. So they should
           start from random coordinates. Defaults to False.
        :type use_random_coordinates: bool, optional
        """

        # Save some values to the object.
        self.smiles = mol.smiles()
        self.orig_smi = mol.orig_smi

        # Set before any of the failure paths below: callers read .energy
        # without first checking .mol (the minimization and ring-conformer
        # steps both do), and add_conformers sorts on it. Infinity is the same
        # placeholder used when the force field cannot score a conformer, so a
        # failed one loses every comparison.
        self.energy = float("inf")

        # Set for the same reason as energy, and with the values that make a
        # failed conformer inert: minimize() returns immediately rather than
        # handing False to the force field (whose own error handler would then
        # raise on Chem.MolToSmiles(False)), and no atom is offered for
        # alignment.
        self.minimized = True
        self.ids_hvy_atms: list[int] = []

        # A caller can hand us a MyMol whose rdkit_mol is None (e.g. failed
        # reprotonation). deepcopy(None) is None, and the subsequent
        # RemoveAllConformers() would raise AttributeError. Treat it as a
        # failed conformer so add_conformers() skips it.
        if mol.rdkit_mol is None:
            self.mol = False
            self.coord_3d_err_warning(None)
            return

        self.mol = copy.deepcopy(mol.rdkit_mol)

        # Remove any previous conformers.
        self.mol.RemoveAllConformers()

        if conformer is None:
            # The user is providing no conformer. So we must generate it.

            # Note that I have confirmed that the below respects chirality.
            # params is a list of ETKDGv2 parameters generated by this command
            # Description of these parameters can be found at
            # help(AllChem.EmbedMolecule)

            try:
                # Newest version
                # print("HERE")
                params = AllChem.ETKDGv3()
            except Exception:
                try:
                    # Try to use ETKDGv2, but it is only present in the python 3.6
                    # version of RDKit.
                    params = AllChem.ETKDGv2()
                except Exception:
                    # Use the original version of ETKDG if python 2.7 RDKit. This
                    # may be resolved in next RDKit update so we encased this in a
                    # try statement.
                    params = AllChem.ETKDG()

            # The default, but just a sanity check.
            params.enforceChirality = True

            # Set a max number of times it will try to calculate the 3D
            # coordinates. Will save a little time. This should be the default
            # (0) but lets set it anyway
            params.maxIterations = 0

            # Also set whether to start from random coordinates.
            params.useRandomCoords = use_random_coordinates

            # RDKit's embedding draws from its own generator, which
            # random.seed() cannot reach; its default (-1) means a fresh,
            # unrecoverable seed on every call. Derive the seed from the
            # (optionally seeded) Python generator instead, so --random_seed
            # also fixes the coordinates.
            params.randomSeed = random.randint(1, 2**31 - 1)

            # AllChem.EmbedMolecule uses geometry to create inital molecule
            # coordinates. This sometimes takes a very long time.
            try:
                AllChem.EmbedMolecule(self.mol, params)
            except RuntimeError as e:
                self.mol = False
                self.coord_3d_err_warning(e)

            # On rare occasions, the new conformer generating algorithm fails
            # because params.useRandomCoords = False. So if it fails, try
            # again with True.
            if (
                self.mol is not False
                and self.mol.GetNumConformers() == 0
                and use_random_coordinates == False
            ):
                params.useRandomCoords = True
                try:
                    AllChem.EmbedMolecule(self.mol, params)
                except RuntimeError as e:
                    self.mol = False
                    self.coord_3d_err_warning(e)

            # On very rare occasions, the new conformer generating algorithm
            # fails. For example, COC(=O)c1cc(C)nc2c(C)cc3[nH]c4ccccc4c3c12 .
            # In this case, the old one still works. So if no coordinates are
            # assigned, try that one. Parameters must have second_embed set to
            # True for this to happen.
            if (
                self.mol is not False
                and second_embed == True
                and self.mol.GetNumConformers() == 0
            ):
                try:
                    # This legacy call takes no EmbedParameters, so it needs
                    # the seed passed separately or it falls back to RDKit's
                    # unseeded default. It reuses the seed drawn above, so it
                    # must differ from the earlier attempts in some other way
                    # or it would fail exactly as they did: random starting
                    # coordinates (the non-random start is what attempt 2
                    # rescues) and tolerance of triangle-smoothing failures,
                    # a common cause of the failures that remain.
                    AllChem.EmbedMolecule(
                        self.mol,
                        useRandomCoords=True,
                        randomSeed=params.randomSeed,
                        ignoreSmoothingFailures=True,
                    )
                except RuntimeError as e:
                    self.mol = False
                    self.coord_3d_err_warning(e)

            # On rare occasions, both methods fail. For example,
            # O=c1cccc2[C@H]3C[NH2+]C[C@@H](C3)Cn21 Another example:
            # COc1cccc2c1[C@H](CO)[N@H+]1[C@@H](C#N)[C@@H]3C[C@@H](C(=O)[O-])[C@H]([C@H]1C2)[N@H+]3C
            if self.mol is not False and self.mol.GetNumConformers() == 0:
                self.mol = False
                self.coord_3d_err_warning(None)
        else:
            # The user has provided a conformer. Just add it.
            conformer.SetId(0)
            self.mol.AddConformer(conformer, assignId=True)

        # Calculate some energies, other housekeeping.
        if self.mol is not False:
            try:
                ff = AllChem.UFFGetMoleculeForceField(self.mol)
                self.energy = ff.CalcEnergy()
            except Exception:
                utils.log(
                    "Warning: Could not calculate energy for molecule "
                    + Chem.MolToSmiles(self.mol)
                )
                # Example of smiles that cause problem here without try...catch:
                # NC1=NC2=C(N[C@@H]3[C@H](N2)O[C@@H](COP(O)(O)=O)C2=C3S[Mo](S)(=O)(=O)S2)C(=O)N1
                # Every consumer of this attribute sorts on it, and a strained
                # or large ligand can exceed any finite placeholder, so a
                # conformer with no energy at all would then outrank a real
                # one. Infinity sorts last unconditionally.
                self.energy = float("inf")
            self.minimized = False
            self.ids_hvy_atms = [
                a.GetIdx() for a in self.mol.GetAtoms() if a.GetAtomicNum() != 1
            ]

    def coord_3d_err_warning(self, err):
        utils.log(
            f'WARNING: RDKit failed to generate 3D coordinates for a molecule originating from "{self.orig_smi}". The SMILES string of the problematic variant is "{self.smiles}". The variant will be skipped. Specific RDKit error: {err}'
        )

    def conformer(self, conf=None):
        """Get or set the conformer. An optional variable can specify the
           conformer to set. If not specified, this function acts as a get for
           the conformer.

        :param conf: The conformer to set, defaults to None
        :param conf: rdkit.Conformer, optional
        :return: An rdkit.Conformer object, if conf is not specified.
        :rtype: rdkit.Conformer
        """

        if conf is None:
            return self.mol.GetConformers()[0]

        self.mol.RemoveAllConformers()
        self.mol.AddConformer(conf)

    def minimize(self):
        """Minimize (optimize) the geometry of the current conformer if it
        hasn't already been optimized."""

        if self.minimized == True:
            # Already minimized. Don't do it again.
            return

        # Perform the minimization, and save the energy.
        try:
            ff = AllChem.UFFGetMoleculeForceField(self.mol)
            ff.Minimize()
            self.energy = ff.CalcEnergy()
        except Exception:
            utils.log(
                "Warning: Could not calculate energy for molecule "
                + Chem.MolToSmiles(self.mol)
            )
            # Same reasoning as in the constructor: an unknown energy has to
            # lose every comparison, which no finite value can guarantee.
            self.energy = float("inf")
        self.minimized = True

    def align_to_me(self, other_conf):
        """Align another conformer to this one.

        :param other_conf: The other conformer to align.
        :type other_conf: MyConformer
        :return: The aligned MyConformer object.
        :rtype: MyConformer
        """

        # Add the conformer of the other MyConformer object.
        self.mol.AddConformer(other_conf.conformer(), assignId=True)

        # Align them.
        AllChem.AlignMolConformers(self.mol, atomIds=self.ids_hvy_atms)

        # Reset the conformer of the other MyConformer object.
        last_conf = self.mol.GetConformers()[-1]
        other_conf.conformer(last_conf)

        # Remove the added conformer.
        self.mol.RemoveConformer(last_conf.GetId())

        # Return that other object.
        return other_conf

    def MolToMolBlock(self):
        """Prints out the first 500 letters of the molblock version of this
        conformer. Good for debugging."""

        mol_copy = copy.deepcopy(self.mol)  # Use it as a template.
        mol_copy.RemoveAllConformers()
        mol_copy.AddConformer(self.conformer())
        utils.log(Chem.MolToMolBlock(mol_copy)[:500])

    def rmsd_to_me(self, other_conf: "MyConformer") -> float:
        """Calculate the rms distance between this conformer and another one.

        :param other_conf: The other conformer to align.
        :type other_conf: MyConformer
        :return: The rmsd, a float.
        :rtype: float
        """

        # Compare inside the molecule the coordinates actually belong to. A
        # molecule rebuilt from the canonical SMILES is a differently ordered
        # graph, so atom i of it would receive the position of some other
        # source atom, and the deprotonation below would then drop whichever
        # atoms are hydrogens in that rebuilt ordering rather than the
        # hydrogens of the coordinate data.
        if self.mol is False or other_conf.mol is False:
            return float("inf")

        probe = copy.deepcopy(self.mol)
        probe.RemoveAllConformers()
        probe.AddConformer(self.conformer(), assignId=True)
        probe.AddConformer(other_conf.conformer(), assignId=True)

        # Deprotonate so the comparison is over heavy atoms only. This can
        # fail and return None; infinity is the conservative answer for an
        # RMSD that could not be computed, because the caller deduplicates on
        # rmsd <= cutoff, so both conformers are kept rather than one being
        # silently discarded.
        probe = MOH.try_deprotanation(probe)
        if probe is None:
            return float("inf")
        return AllChem.GetConformerRMS(probe, 0, 1, prealigned=True)

    def coords(self):
        """Get the coordinates of this conformer. For debugging.

        :return: A list of coordinates.
        :rtype: list
        """

        return self.conformer().GetPositions()

    def write_pdb_file(self, filename: str):
        """Write this conformer to a PDB file. For debugging.

        :param filename: The name of the file to write.
        :type filename: str
        """

        # Make a new molecule.
        mol = copy.deepcopy(self.mol)
        mol.RemoveAllConformers()

        # Add the conformer of the other MyConformer object.
        mol.AddConformer(self.conformer(), assignId=True)

        # Write the PDB file.
        AllChem.MolToPDBFile(mol, filename)

    def get_energy(self) -> float:
        """Get the energy of this conformer. For debugging.

        :return: The energy.
        :rtype: float
        """

        ff = AllChem.UFFGetMoleculeForceField(self.mol)
        return ff.CalcEnergy()
