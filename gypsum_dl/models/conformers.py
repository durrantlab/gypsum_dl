import copy
import operator

from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem

from gypsum_dl import handlers


class Conformer:
    """
    Encapsulates a single RDKit Conformer object and its associated energy.
    """

    def __init__(
        self,
        mol: Chem.Mol,
        rdkit_conformer: Chem.Conformer,
        energy: float,
        minimized: bool = False,
    ):
        self._mol = mol  # Keep a reference to the parent RDKit Mol for atom access
        self._rdkit_conformer = rdkit_conformer
        self._energy = energy
        self._minimized = minimized
        self._heavy_atom_ids = [
            a.GetIdx() for a in self._mol.GetAtoms() if a.GetAtomicNum() != 1
        ]

    @property
    def rdkit_conformer(self) -> Chem.Conformer:
        return self._rdkit_conformer

    @property
    def energy(self) -> float:
        return self._energy

    @property
    def is_minimized(self) -> bool:
        return self._minimized

    def minimize(self) -> None:
        """
        Minimizes the geometry of the current conformer if it hasn't already been optimized.
        """
        if self._minimized:
            return

        try:
            # Create a temporary Mol just for minimization of this conformer
            temp_mol = Chem.Mol(self._mol)
            temp_mol.RemoveAllConformers()
            temp_mol.AddConformer(self._rdkit_conformer, assignId=True)

            ff = AllChem.UFFGetMoleculeForceField(temp_mol)  # type: ignore
            ff.Minimize()
            self._energy = ff.CalcEnergy()
            self._rdkit_conformer = temp_mol.GetConformers()[
                0
            ]  # Update with minimized conformer
            self._minimized = True
            logger.info(f"Conformer minimized to energy: {self._energy:.2f}")
        except Exception as e:
            logger.warning(f"Could not minimize conformer. Error: {e}")
            self._energy = 9999.0  # Assign a high energy to mark as problematic

    def align_to_me(self, other_conformer: "Conformer") -> None:
        """
        Aligns another conformer to this conformer in place.
        """
        # Create a temporary mol with both conformers to align
        temp_mol = Chem.Mol(self._mol)
        temp_mol.RemoveAllConformers()
        temp_mol.AddConformer(self.rdkit_conformer, assignId=True)
        temp_mol.AddConformer(other_conformer.rdkit_conformer, assignId=True)

        AllChem.AlignMolConformers(temp_mol, atomIds=self._heavy_atom_ids)  # type: ignore

        # Update the other conformer's RDKit Conformer object
        other_conformer._rdkit_conformer = temp_mol.GetConformers()[1]

    def rmsd_to_me(self, other_conformer: "Conformer") -> float:
        """
        Calculates the RMSD between this conformer and another one.
        Assumes both conformers belong to the same parent molecule structure.
        """
        # A new mol is created to ensure atom indices match and for RMSD calculation
        temp_mol = Chem.Mol(self._mol)
        temp_mol.RemoveAllConformers()
        temp_mol.AddConformer(self.rdkit_conformer, assignId=True)
        temp_mol.AddConformer(other_conformer.rdkit_conformer, assignId=True)

        # The handlers.try_deprotanation might be problematic if it changes atom order.
        # For RMSD calculation, it's generally best to ensure consistent atom indexing
        # between the two conformers you're comparing. If the mol object itself is modified,
        # it can break the RMSD calculation. Assuming the mol object passed to Conformer
        # and used by rmsd_to_me always represents the same atom order for the comparison.
        mol_for_rmsd = handlers.try_deprotanation(temp_mol)
        if mol_for_rmsd is None:
            raise RuntimeError(
                "Failed to prepare molecule for RMSD calculation (deprotonation failed)."
            )
        return AllChem.GetConformerRMS(mol_for_rmsd, 0, 1, prealigned=True)

    def get_positions(self) -> list[list[float]]:
        """Returns the 3D coordinates of the conformer."""
        return self._rdkit_conformer.GetPositions().tolist()


class ConformerSet:
    """
    Manages a collection of conformers for a specific molecule.
    """

    def __init__(self, molecule: Chem.Mol):
        self._molecule_template = (
            molecule  # Store a template Mol for generating new conformers
        )
        self._conformers: list[Conformer] = []

    def generate_conformers(
        self,
        num_to_generate: int,
        use_random_coordinates: bool = False,
        embed_second_try: bool = False,
    ) -> None:
        """
        Generates new conformers for the associated molecule.
        """
        initial_count = len(self._conformers)
        num_needed = num_to_generate - initial_count

        for _ in range(num_needed):
            mol_copy = Chem.Mol(self._molecule_template)  # Work on a copy
            mol_copy.RemoveAllConformers()  # Ensure no existing conformers

            params = None
            try:
                params = AllChem.ETKDGv3()  # type: ignore
            except AttributeError:  # Fallback for older RDKit versions
                try:
                    params = AllChem.ETKDGv2()  # type: ignore
                except AttributeError:
                    params = AllChem.ETKDG()  # type: ignore

            params.enforceChirality = True
            params.maxIterations = 0
            params.useRandomCoords = use_random_coordinates

            try:
                AllChem.EmbedMolecule(mol_copy, params)  # type: ignore
            except Exception as e:
                logger.warning(
                    f"Initial conformer embedding failed for {Chem.MolToSmiles(mol_copy)}. Error: {e}"
                )
                if mol_copy.GetNumConformers() == 0 and embed_second_try:
                    try:
                        AllChem.EmbedMolecule(  # type: ignore
                            mol_copy, useRandomCoords=use_random_coordinates
                        )
                    except Exception as e_second:
                        logger.warning(
                            f"Second conformer embedding failed. Error: {e_second}"
                        )

            if mol_copy.GetNumConformers() > 0:
                rdkit_conf = mol_copy.GetConformers()[0]
                try:
                    ff = AllChem.UFFGetMoleculeForceField(mol_copy)  # type: ignore
                    energy = ff.CalcEnergy()
                except Exception as e:
                    logger.warning(
                        f"Could not calculate energy for generated conformer. Error: {e}"
                    )
                    energy = 9999.0
                self._conformers.append(
                    Conformer(self._molecule_template, rdkit_conf, energy)
                )
            else:
                logger.warning("Failed to generate a conformer.")

        self.sort_conformers()

    def minimize_all(self) -> None:
        """Minimizes all conformers in the set."""
        for conf in self._conformers:
            conf.minimize()
        self.sort_conformers()

    def eliminate_structurally_similar(self, rmsd_cutoff: float = 0.1) -> None:
        """
        Eliminates conformers that are very geometrically similar based on RMSD.
        Assumes conformers are already sorted by energy for better selection.
        """
        if not self._conformers:
            return

        unique_conformers: list[Conformer] = []
        if self._conformers:
            unique_conformers.append(self._conformers[0])

        for i in range(1, len(self._conformers)):
            is_unique = True
            current_conf = self._conformers[i]
            for unique_conf in unique_conformers:
                # Align current_conf to unique_conf before calculating RMSD
                aligned_current_conf = copy.deepcopy(
                    current_conf
                )  # Avoid modifying original during check
                unique_conf.align_to_me(aligned_current_conf)
                rmsd = unique_conf.rmsd_to_me(aligned_current_conf)
                if rmsd <= rmsd_cutoff:
                    is_unique = False
                    break
            if is_unique:
                unique_conformers.append(current_conf)

        self._conformers = unique_conformers
        logger.info(
            f"Reduced to {len(self._conformers)} unique conformers after RMSD cutoff of {rmsd_cutoff}."
        )
        self.sort_conformers()

    def sort_conformers(self) -> None:
        """Sorts conformers by energy."""
        self._conformers.sort(key=operator.attrgetter("energy"))

    def get_conformers(self) -> list[Conformer]:
        """Returns the list of managed Conformer objects."""
        return self._conformers

    def load_into_rdkit_mol(self, target_mol: Chem.Mol) -> None:
        """Loads all stored conformers into a given RDKit Mol object."""
        target_mol.RemoveAllConformers()
        for conf_obj in self._conformers:
            target_mol.AddConformer(conf_obj.rdkit_conformer, assignId=True)
        logger.info(f"Loaded {len(self._conformers)} conformers into RDKit Mol object.")
