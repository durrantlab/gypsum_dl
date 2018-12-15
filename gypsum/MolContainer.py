import __future__


import gypsum.Utils as Utils
import gypsum.ChemUtils as ChemUtils
import gypsum.MyMol as MyMol

try:
    from molvs import standardize_smiles as ssmiles
except:
    Utils.log("You need to install molvs and its dependencies.")
    raise ImportError("You need to install molvs and its dependencies.")

try:
    from rdkit import Chem
except:
    Utils.log("You need to install rdkit and its dependencies.")
    raise ImportError("You need to install rdkit and its dependencies.")

class MolContainer:
    """
    The molecucle container class. It stores all the molecules (tautomers,
    etc.) associated with a single input SMILES entry.
    """

    def __init__(self, smiles, name, index, properties):
        """The constructor for this class.

        :param [str] smiles: A list of SMILES strings.

        :param [str] name: The name of the molecule.

        :param int index: The index of this MolContainer in the main
                          MolContainer list.
        :param {dict} properties: A dictionary of properties from the sdf.
        """

        # Set some variables on the contnr level (not the MyMolecule level)
        self.contnr_idx = index
        self.orig_smi = smiles
        self.orig_smi_deslt = smiles  # initial assumption
        self.mols = []
        #self.smiles_to_pH = {}
        self.name = name
        self.properties = properties

        self.mol_orig_frm_inp_smi = MyMol.MyMol(smiles, name)
        self.mol_orig_frm_inp_smi.contnr_idx = self.contnr_idx

        self.frgs = ""  # For caching.

        # Save the original canonical smiles
        self.orig_smi_canonical = self.mol_orig_frm_inp_smi.smiles()

        # Get the number of nonaromatic rings
        self.num_nonaro_rngs = len(self.mol_orig_frm_inp_smi.get_idxs_of_nonaro_rng_atms())

        # Get the number of chiral centers, assigned or not
        self.num_specif_chiral_cntrs = len(
            self.mol_orig_frm_inp_smi.chiral_cntrs_only_asignd()
        )

        # Get the non-acidic carbon-hydrogen footprint.
        self.carbon_hydrogen_count = self.mol_orig_frm_inp_smi.count_hyd_bnd_to_carb()

    def mol_with_smiles_is_in_container(self, smiles):
        """
        Checks whether or not a given smiles string is already in this
        container.

        :param str smiles: The smiles string to check.

        :returns: True if it is present, otherwise a new MyMol.MyMol object
                  corresponding to that smiles.
        :rtype: :class:`str` ???
        """

        # Checks all the mols in this container to see if a given smiles is
        # already present. Returns a new MyMol object if it isn't, True
        # otherwise.

        # First, get the set of all cannonical smiles
        # TODO: Probably shouldn't be generating this on the fly every time
        # you use it!
        can_smi_in_this_container = set([m.smiles() for m in self.mols])

        # Get this new smiles in cannonical form
        amol = MyMol.MyMol(smiles)
        if amol.smiles() in can_smi_in_this_container:
            return True
        else:
            return amol

    def add_smiles(self, smiles):
        """
        Adds smiles strings to this container. SMILES are always isomeric and
        always unique (canonical).

        :param str|[str] smiles: A list of SMILES strings. If it's a string,
                         it is converted into a list.
        """

        if isinstance(smiles, str):  # smiles must be array of strs
            smiles = [smiles]

        # Keep only the mols with smiles that are not already present.
        for s in smiles:
            result = self.mol_with_smiles_is_in_container(s)
            if result != True:
                # Much of the contnr info should be passed to each
                # molecule, too, for convenience.
                result.name = self.name
                result.name = self.orig_smi
                result.orig_smi_canonical = self.orig_smi_canonical
                result.orig_smi_deslt = self.orig_smi_deslt
                result.contnr_idx = self.contnr_idx

                self.mols.append(result)

    def add_mol(self, mol):
        """Adds a molecule to this container. Does NOT check for uniqueness.

        :param MyMol.MyMol mol: The MyMol.MyMol object to add.
        """

        self.mols.append(mol)

    def all_smiles(self):
        """
        Gets a list of all the noh canonical smiles in this container.

        :returns: a [str], the canonical, noh smiles string.
        :rtype: :class:`str` ???
        """
        smiles = []
        for m in self.mols:
            if m.rdkit_mol is not None:
                smiles.append(m.smiles())  # Pass true as a parameter to skip hydrogens.

        return smiles

    def get_frags_of_orig_smi(self):
        """
        Gets a list of the fragments found in the original smiles string
        passed to this container.

        :returns: a list of the fragments, as rdkit.Mol objects.
        :rtype: :class:`str` ???
        """

        if self.frgs != "":
            return self.frgs

        frags = self.mol_orig_frm_inp_smi.get_frags_of_orig_smi()
        self.frgs = frags
        return frags

    def update_orig_smi(self, orig_smi):
        """
        Updates the orig_smi string.

        :param str orig_smi: The replacement smiles string.
        """

        # Update the MolContainer object
        self.orig_smi = orig_smi
        self.orig_smi_deslt = orig_smi
        self.mol_orig_frm_inp_smi = MyMol.MyMol(self.orig_smi, self.name)
        self.frgs = ""
        self.orig_smi_canonical = self.mol_orig_frm_inp_smi.smiles()
        self.num_nonaro_rngs = len(self.mol_orig_frm_inp_smi.get_idxs_of_nonaro_rng_atms())
        self.num_specif_chiral_cntrs = len(
            self.mol_orig_frm_inp_smi.chiral_cntrs_only_asignd()
        )

        # None of the mols derived to date, if present, are accurate.
        self.mols = []

    def add_container_properties(self):
        """
        Adds all properties from the container, which is currently populated
        from the SDF when using 2D-SDFs as input, to the molecules.
        """
        for mol in self.mols:
            mol.mol_props.update(self.properties)
            mol.set_all_rdkit_mol_props()

    def remove_identical_mols_from_container(self):
        # For reasons I don't understand, the following doesn't give unique
        # canonical smiles:

        # Chem.MolToSmiles(self.mols[0].rdkit_mol, isomericSmiles=True,
        # canonical=True)

        # This block for debugging
        all_smiles = [m.smiles() for m in self.mols]  # Get all the smiles as stored.

        wrong_cannonical_smiles = [
            Chem.MolToSmiles(
                m.rdkit_mol,  # Using the RdKit mol stored in MyMol
                isomericSmiles=True,
                canonical=True
            ) for m in self.mols
        ]

        right_cannonical_smiles = [
            Chem.MolToSmiles(
                Chem.MolFromSmiles(  # Regenerating the RdKit mol from the smiles string stored in MyMol
                    m.smiles()
                ),
                isomericSmiles=True,
                canonical=True
            ) for m in self.mols]

        if len(set(wrong_cannonical_smiles)) != len(set(right_cannonical_smiles)):
            Utils.log("ERROR!")
            Utils.log("Stored smiles string in this container:")
            Utils.log("\n".join(all_smiles))
            Utils.log("")
            Utils.log("""Supposedly cannonical smiles strings generated from stored
                RDKit Mols in this container:""")
            Utils.log("\n".join(wrong_cannonical_smiles))
            Utils.log("""But if you plop these into chemdraw, you'll see some of them
                represent identical structures.""")
            Utils.log("")
            Utils.log("""Cannonical smiles strings generated from RDKit mols that
                were generated from the stored smiles string in this container:""")
            Utils.log("\n".join(right_cannonical_smiles))
            Utils.log("""Now you see the identical molecules. But why didn't the previous
                method catch them?""")
            Utils.log("")

            Utils.log("""Note that the third method identifies duplicates that the second
                method doesn't.""")
            Utils.log("")
            Utils.log("=" * 20)

        # You need to make new molecules to get it to work.
        new_smiles = [m.smiles() for m in self.mols]
        new_mols = [Chem.MolFromSmiles(smi) for smi in new_smiles]
        new_can_smiles = [Chem.MolToSmiles(new_mol, isomericSmiles=True, canonical=True) for new_mol in new_mols]

        can_smiles_already_set = set([])
        for i, new_can_smile in enumerate(new_can_smiles):
            if not new_can_smile in can_smiles_already_set:
                # Never seen before
                can_smiles_already_set.add(new_can_smile)
            else:
                # See before. Delete!
                self.mols[i] = None

        while None in self.mols:
            self.mols.remove(None)

    def update_idx(self, new_idx):
        if type(new_idx)!= int:
            raise Exception("new idx value must be an int")
        self.contnr_idx = new_idx
        self.mol_orig_frm_inp_smi.contnr_idx = self.contnr_idx
