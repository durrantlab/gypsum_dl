"""
This module describes the MoleculeContainer, which contains different Molecule
objects. Each object in this container is derived from the same input molecule
(so they are variants). Note that conformers (3D coordinate sets) live inside
Molecule.
"""

from typing import Any

from gypsum_dl import Molecule, chem_utils


class MoleculeContainer:
    """The molecule container class. It stores all the molecules (tautomers,
    etc.) associated with a single input SMILES entry."""

    def __init__(
        self, smiles: str, name: str, index: int, properties: dict[str, Any]
    ) -> None:
        """The constructor.

        Args:
            smiles: A list of SMILES strings.
            name: The name of the molecule.
            index: The index of this MoleculeContainer in the main MoleculeContainer
                list.
            properties: A dictionary of properties from the sdf.
        """

        # Set some variables are set on the container level (not the Molecule
        # level)
        self.contnr_idx = index
        self.contnr_idx_orig = index  # Because if some circumstances (mpi),
        # might be reset. But good to have
        # original for filename output.
        self.orig_smi = smiles
        self.orig_smi_deslt = smiles  # initial assumption
        self.mols = []
        self.name = name
        self.properties = properties
        self.mol_orig_frm_inp_smi = Molecule(smiles, name)
        self.mol_orig_frm_inp_smi.contnr_idx = self.contnr_idx
        self.frgs = ""  # For caching.

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

        # Also get the number of chiral centers, unassigned
        self.num_unspecif_chiral_cntrs = len(
            self.mol_orig_frm_inp_smi.chiral_cntrs_w_unasignd()
        )

        # Get the non-acidic carbon-hydrogen footprint.
        self.carbon_hydrogen_count = self.mol_orig_frm_inp_smi.count_hyd_bnd_to_carb()

    def mol_with_smiles_is_in_contnr(self, smiles: str) -> bool | Molecule:
        """Checks whether or not a given smiles string is already in this
           container.

        Args:
            smiles: The smiles string to check.

        Returns:
            True if it is present, otherwise a new Molecule object
                corresponding to that smiles.
        """

        # Checks all the mols in this container to see if a given smiles is
        # already present. Returns a new Molecule object if it isn't, True
        # otherwise.

        # First, get the set of all cannonical smiles.
        # TODO: Probably shouldn't be generating this on the fly every time
        # you use it!
        can_smi_in_this_container = {m.smiles() for m in self.mols}

        # Determine whether it is already in the container, and act
        # accordingly.
        amol = Molecule(smiles)
        return True if amol.smiles() in can_smi_in_this_container else amol

    def add_smiles(self, smiles: list[str]) -> None:
        """Adds smiles strings to this container. SMILES are always isomeric
           and always unique (canonical).

        Args:
            smiles: A list of SMILES strings. If it's a string, it is
                converted into a list.
        """

        # Convert it into a list if it comes in as a string.
        if isinstance(smiles, str):
            smiles = [smiles]

        # Keep only the mols with smiles that are not already present.
        for s in smiles:
            result = self.mol_with_smiles_is_in_contnr(s)
            if not result:
                # Much of the contnr info should be passed to each molecule,
                # too, for convenience.
                result.name = self.name
                result.name = self.orig_smi
                result.orig_smi_canonical = self.orig_smi_canonical
                result.orig_smi_deslt = self.orig_smi_deslt
                result.contnr_idx = self.contnr_idx

                self.mols.append(result)

    def add_mol(self, mol: Molecule) -> None:
        """Adds a molecule to this container. Does NOT check for uniqueness.

        Args:
            mol: The Molecule object to add.
        """

        self.mols.append(mol)

    def all_can_noh_smiles(self) -> str:
        """Gets a list of all the noh canonical smiles in this container.

        Args:
            The canonical, noh smiles string.
        """

        # True means noh
        return [m.smiles(True) for m in self.mols if m.rdkit_mol is not None]

    def get_frags_of_orig_smi(self) -> list:
        """Gets a list of the fragments found in the original smiles string
           passed to this container.

        Returns:
            A list of the fragments, as rdkit.Mol objects. Also saves to self.frgs.
        """

        if self.frgs != "":
            return self.frgs

        frags = self.mol_orig_frm_inp_smi.get_frags_of_orig_smi()
        self.frgs = frags
        return frags

    def update_orig_smi(self, orig_smi: str) -> None:
        """Updates the orig_smi string. Used by desalter (to replace with
           largest fragment).

        Args:
            orig_smi: The replacement smiles string.
        """

        # Update the MoleculeContainer object
        self.orig_smi = orig_smi
        self.orig_smi_deslt = orig_smi
        self.mol_orig_frm_inp_smi = Molecule(self.orig_smi, self.name)
        self.frgs = ""
        self.orig_smi_canonical = self.mol_orig_frm_inp_smi.smiles()
        self.num_nonaro_rngs = len(
            self.mol_orig_frm_inp_smi.get_idxs_of_nonaro_rng_atms()
        )
        self.num_specif_chiral_cntrs = len(
            self.mol_orig_frm_inp_smi.chiral_cntrs_only_asignd()
        )
        self.num_unspecif_chiral_cntrs = len(
            self.mol_orig_frm_inp_smi.chiral_cntrs_w_unasignd()
        )

        # None of the mols derived to date, if present, are accurate.
        self.mols = []

    def add_container_properties(self) -> None:
        """Adds all properties from the container to the molecules. Used when
        saving final files, to keep a record in the file itself."""

        for mol in self.mols:
            mol.mol_props.update(self.properties)
            mol.set_all_rdkit_mol_props()

    def remove_identical_mols_from_contnr(self) -> None:
        """Removes identical molecules from this container."""

        # For reasons I don't understand, the following doesn't give unique
        # canonical smiles:

        # Chem.MolToSmiles(self.mols[0].rdkit_mol, isomericSmiles=True,
        # canonical=True)

        # # This block for debugging. JDD: Needs attention?
        # all_can_noh_smiles = [m.smiles() for m in self.mols]  # Get all the smiles as stored.

        # wrong_cannonical_smiles = [
        #     Chem.MolToSmiles(
        #         m.rdkit_mol,  # Using the RdKit mol stored in Molecule
        #         isomericSmiles=True,
        #         canonical=True
        #     ) for m in self.mols
        # ]

        # right_cannonical_smiles = [
        #     Chem.MolToSmiles(
        #         Chem.MolFromSmiles(  # Regenerating the RdKit mol from the smiles string stored in Molecule
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

    def update_idx(self, new_idx: int) -> None:
        """Updates the index of this container.

        Args:
            new_idx: The new index.
        """

        if not isinstance(new_idx, int):
            raise TypeError("New idx value must be an int.")
        self.contnr_idx = new_idx
        self.mol_orig_frm_inp_smi.contnr_idx = self.contnr_idx
