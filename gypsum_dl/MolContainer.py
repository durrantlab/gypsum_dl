"""
This module describes the MolContainer, which contains different MyMol.MyMol
objects. Each object in this container is derived from the same input molecule
(so they are variants). Note that conformers (3D coordinate sets) live inside
MyMol.MyMol. So, just to clarify:

MolContainer.MolContainer > MyMol.MyMol > MyMol.MyConformers
"""

import copy

from gypsum_dl import MyMol, chem_utils, utils


class MolContainer:
    """The molecucle container class. It stores all the molecules (tautomers,
    etc.) associated with a single input SMILES entry."""

    def __init__(self, smiles, name, index, properties):
        """The constructor.

        :param smiles: A list of SMILES strings.
        :type smiles: str
        :param name: The name of the molecule.
        :type name: str
        :param index: The index of this MolContainer in the main MolContainer
           list.
        :type index: int
        :param properties: A dictionary of properties from the sdf.
        :type properties: dict
        """

        # Set some variables are set on the container level (not the MyMol
        # level)
        self.contnr_idx = index
        self.contnr_idx_orig = index  # Because if some circumstances (mpi),
        # might be reset. But good to have
        # original for filename output.
        self.orig_smi = smiles
        # The string the user submitted. orig_smi and orig_smi_deslt are both
        # overwritten with the largest fragment by update_orig_smi, so neither
        # can serve as the record of the input once a salt has been stripped.
        self.orig_smi_input = smiles
        self.orig_smi_deslt = smiles  # initial assumption
        self.mols = []
        self.name = name
        self.properties = properties

        # Ranking energies for this compound's variants, keyed on canonical
        # SMILES. Every SMILES step prunes its variants by embedding a
        # throwaway conformer per candidate, and a variant that survives one
        # step is a candidate again in the next, so without a cache the same
        # structure is embedded once per step. Scoped to the container rather
        # than the module because that is where the repeats are (all variants
        # of one compound) and because it bounds the cache to a compound's
        # variant count instead of the whole input library. See
        # chem_utils.first_conf_energy.
        self.probe_energies: dict[str, float | None] = {}

        # Everything derived from orig_smi (the reference molecule, its
        # canonical smiles, the ring/chiral counts, and the fragment cache) is
        # built in one place.
        self.derive_from_orig_smi()

    def derive_from_orig_smi(self) -> None:
        """Rebuild every field that is a function of self.orig_smi.

        The constructor and update_orig_smi both need this derivation, and
        keeping two copies of it let them drift: one path once refreshed a
        derived field that the other left describing the pre-desalt molecule.
        One copy means a derived field cannot be added to one path and
        forgotten on the other.

        Returns:
            None. Sets the derived attributes on this container.
        """

        self.mol_orig_frm_inp_smi = MyMol.MyMol(self.orig_smi, self.name)
        self.mol_orig_frm_inp_smi.contnr_idx = self.contnr_idx
        self.frgs = MyMol.UNSET  # For caching.

        # Save the original canonical smiles
        self.orig_smi_canonical = self.mol_orig_frm_inp_smi.smiles()

        # Get the number of nonaromatic rings
        self.num_nonaro_rngs = len(
            self.mol_orig_frm_inp_smi.get_idxs_of_nonaro_rng_atms()
        )

        # Get the number of chiral centers, assigned
        self.num_specif_chiral_cntrs = len(
            self.mol_orig_frm_inp_smi.chiral_cntrs_only_asignd()
        )

        # Get the total number of chiral centers, assigned or not. The name
        # says unspecified, but chiral_cntrs_w_unasignd returns both kinds, so
        # this count includes the assigned centers tallied above.
        self.num_unspecif_chiral_cntrs = len(
            self.mol_orig_frm_inp_smi.chiral_cntrs_w_unasignd()
        )

    def contnr_props(self) -> "MyMol.ContnrProps":
        """Collect the container-level fields every variant here shares.

        Lets a step hand a job the few fields that describe the input compound
        instead of this container, which holds every variant of it. See
        MyMol.ContnrProps.

        Returns:
            This container's copy of the container-level fields.
        """

        return {
            "contnr_idx": self.contnr_idx,
            "name": self.name,
            "orig_smi": self.orig_smi,
            "orig_smi_deslt": self.orig_smi_deslt,
            "orig_smi_canonical": self.orig_smi_canonical,
        }

    def copy_of_orig_mol(self) -> "MyMol.MyMol":
        """Hand back an independent copy of this container's reference molecule.

        mol_orig_frm_inp_smi is the container's record of the input structure,
        and several steps fall back to "just use the molecule as given" (the
        single-fragment desalting path, the desalting and ionization failure
        paths). Adding that object itself to self.mols aliases the record into
        the working set, where later steps mutate it in place:
        make_first_3d_conf_no_min replaces its rdkit_mol with a reprotonated,
        3D one. Whether that happens depends on the job manager, since the
        multiprocessing and mpi paths hand the steps pickled copies. A copy
        makes every path behave like the pickling ones and keeps the record
        describing the input.

        Returns:
            A deep copy of the reference molecule, carrying this container's
                index.
        """

        mol_copy = copy.deepcopy(self.mol_orig_frm_inp_smi)
        mol_copy.contnr_idx = self.contnr_idx
        return mol_copy

    def contains_canonical_smiles(self, can_smi: str | None) -> bool:
        """Report whether an already-canonicalized smiles is in this container.

        Split out from mol_with_smiles_is_in_contnr so add_smiles can ask the
        question about a molecule it has already built, instead of building a
        second one to ask with. Non-string arguments (smiles() reports failure
        as None) are never considered present: two molecules that both failed
        to canonicalize have not been shown to be the same molecule.

        Args:
            can_smi: The canonical smiles string to look for.

        Returns:
            True if a molecule with that canonical smiles is already here.
        """

        if not isinstance(can_smi, str):
            return False

        # TODO: Probably shouldn't be generating this on the fly every time
        # you use it!
        return can_smi in {m.smiles() for m in self.mols}

    def mol_with_smiles_is_in_contnr(self, smiles: str) -> bool:
        """Checks whether or not a given smiles string is already in this
           container.

        :param smiles: The smiles string to check.
        :type smiles: str
        :return: True if it is present, otherwise False.
        :rtype: bool
        """

        return self.contains_canonical_smiles(MyMol.MyMol(smiles).smiles())

    def add_smiles(self, smiles):
        """Adds smiles strings to this container. SMILES are always isomeric
           and always unique (canonical).

        :param smiles: A list of SMILES strings. If it's a string, it is
           converted into a list.
        :type smiles: str
        """

        # Convert it into a list if it comes in as a string.
        if isinstance(smiles, str):
            smiles = [smiles]

        # Keep only the mols with smiles that are not already present.
        for s in smiles:
            amol = MyMol.MyMol(s)
            if self.contains_canonical_smiles(amol.smiles()):
                continue

            # Much of the contnr info should be passed to each molecule,
            # too, for convenience. Go through inherit_contnr_props rather than
            # assigning the fields here: this path and that one used to carry
            # separate copies of the list, so a molecule's fields depended on
            # which one had built it.
            amol.inherit_contnr_props(self)

            self.mols.append(amol)

    def add_mol(self, mol):
        """Adds a molecule to this container. Does NOT check for uniqueness.

        :param mol: The MyMol.MyMol object to add.
        :type mol: MyMol.MyMol
        """

        self.mols.append(mol)

    def all_can_noh_smiles(self):
        """Gets a list of all the noh canonical smiles in this container.

        :return: The canonical, noh smiles string.
        :rtype: str
        """

        # True means noh
        return [m.smiles(True) for m in self.mols if m.rdkit_mol is not None]

    def get_frags_of_orig_smi(self):
        """Gets a list of the fragments found in the original smiles string
           passed to this container.

        :return: A list of the fragments, as rdkit.Mol objects. Also saves to
           self.frgs.
        :rtype: list
        """

        if self.frgs is not MyMol.UNSET:
            return self.frgs

        frags = self.mol_orig_frm_inp_smi.get_frags_of_orig_smi()
        self.frgs = frags
        return frags

    def update_orig_smi(self, orig_smi):
        """Updates the orig_smi string. Used by desalter (to replace with
           largest fragment).

        :param orig_smi: The replacement smiles string.
        :type orig_smi: str
        """

        # Update the MolContainer object. orig_smi_input is deliberately left
        # alone: it is what the user submitted, and reports that name the input
        # (the PDB header, the failure file) read from it.
        self.orig_smi = orig_smi
        self.orig_smi_deslt = orig_smi

        # Refresh everything that describes orig_smi; otherwise the derived
        # counts keep describing the pre-desalt (salted) molecule.
        self.derive_from_orig_smi()

        # None of the mols derived to date, if present, are accurate.
        self.mols = []

    def add_container_properties(self):
        """Adds all properties from the container to the molecules. Used when
        saving final files, to keep a record in the file itself."""

        # Input-file properties must not overwrite values Gypsum-DL computed
        # itself (Energy, UniqueID, and so on), so they only fill gaps.
        # Energy needs more than that: when no step computes one (2D output,
        # or both optimization and ring conformations skipped), filling the
        # gap would report the input structure's energy as this variant's.
        # The input value is kept, but under its own tag. It replaces any
        # Input_Energy the input already carried, since the input's Energy
        # describes the structure actually submitted.
        input_props = dict(self.properties)
        if "Energy" in input_props:
            input_props["Input_Energy"] = input_props.pop("Energy")

        for mol in self.mols:
            for key, val in input_props.items():
                mol.mol_props.setdefault(key, val)
            mol.set_all_rdkit_mol_props()

    def remove_identical_mols_from_contnr(self):
        """Removes itentical molecules from this container."""

        # For reasons I don't understand, the following doesn't give unique
        # canonical smiles:

        # Chem.MolToSmiles(self.mols[0].rdkit_mol, isomericSmiles=True,
        # canonical=True)

        # # This block for debugging. JDD: Needs attention?
        # all_can_noh_smiles = [m.smiles() for m in self.mols]  # Get all the smiles as stored.

        # wrong_cannonical_smiles = [
        #     Chem.MolToSmiles(
        #         m.rdkit_mol,  # Using the RdKit mol stored in MyMol
        #         isomericSmiles=True,
        #         canonical=True
        #     ) for m in self.mols
        # ]

        # right_cannonical_smiles = [
        #     Chem.MolToSmiles(
        #         Chem.MolFromSmiles(  # Regenerating the RdKit mol from the smiles string stored in MyMol
        #             m.smiles()
        #         ),
        #         isomericSmiles=True,
        #         canonical=True
        #     ) for m in self.mols]

        # if len(set(wrong_cannonical_smiles)) != len(set(right_cannonical_smiles)):
        #     utils.log("ERROR!")
        #     utils.log("Stored smiles string in this container:")
        #     utils.log("\n".join(all_can_noh_smiles))
        #     utils.log("")
        #     utils.log("""Supposedly cannonical smiles strings generated from stored
        #         RDKit Mols in this container:""")
        #     utils.log("\n".join(wrong_cannonical_smiles))
        #     utils.log("""But if you plop these into chemdraw, you'll see some of them
        #         represent identical structures.""")
        #     utils.log("")
        #     utils.log("""Cannonical smiles strings generated from RDKit mols that
        #         were generated from the stored smiles string in this container:""")
        #     utils.log("\n".join(right_cannonical_smiles))
        #     utils.log("""Now you see the identical molecules. But why didn't the previous
        #         method catch them?""")
        #     utils.log("")

        #     utils.log("""Note that the third method identifies duplicates that the second
        #         method doesn't.""")
        #     utils.log("")
        #     utils.log("=" * 20)

        # # You need to make new molecules to get it to work.
        # new_smiles = [m.smiles() for m in self.mols]
        # new_mols = [Chem.MolFromSmiles(smi) for smi in new_smiles]
        # new_can_smiles = [Chem.MolToSmiles(new_mol, isomericSmiles=True, canonical=True) for new_mol in new_mols]

        # can_smiles_already_set = set([])
        # for i, new_can_smile in enumerate(new_can_smiles):
        #     if not new_can_smile in can_smiles_already_set:
        #         # Never seen before
        #         can_smiles_already_set.add(new_can_smile)
        #     else:
        #         # Seen before. Delete!
        #         self.mols[i] = None

        # while None in self.mols:
        #     self.mols.remove(None)

        self.mols = chem_utils.uniq_mols_in_list(self.mols)

    def update_idx(self, new_idx):
        """Updates the index of this container.

        :param new_idx: The new index.
        :type new_idx: int
        """

        if type(new_idx) != int:
            utils.exception("New idx value must be an int.")
        self.contnr_idx = new_idx
        self.mol_orig_frm_inp_smi.contnr_idx = self.contnr_idx
