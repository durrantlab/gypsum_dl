from loguru import logger

from gypsum_dl.models.conformers import ConformerSet
from gypsum_dl.models.molecule import Molecule


class MoleculeContainer:
    """
    A container to group related Molecule objects (e.g., tautomers, protonation states)
    and manage their preparation workflow, including conformer generation.
    """

    def __init__(self, initial_molecule: Molecule, container_id: int | None = None):
        self._initial_molecule = (
            initial_molecule  # The original input molecule (e.g., desalted)
        )
        self._variants: list[Molecule] = [
            initial_molecule
        ]  # Different tautomers, stereoisomers etc.
        self._container_id = container_id
        # Using dict[str, ConformerSet] for type hinting
        self._conformer_sets: dict[
            str, ConformerSet
        ] = {}  # Key: unique identifier for a variant (e.g., canonical SMILES)

        # Initialize a conformer set for the initial molecule
        initial_mol_key = initial_molecule.canonical_smiles(True)
        self.name: str = initial_molecule.name
        if initial_mol_key is None:
            logger.warning(
                f"Initial molecule {initial_molecule.name} has no canonical SMILES, cannot manage conformers effectively."
            )
        else:
            self._conformer_sets[initial_mol_key] = ConformerSet(
                initial_molecule.rdkit_mol
            )

    @property
    def initial_molecule(self) -> Molecule:
        return self._initial_molecule

    @property
    def container_id(self) -> int | None:
        return self._container_id

    def add_variant(self, molecule_variant: Molecule) -> None:
        """Adds a new molecular variant to this container."""
        if molecule_variant not in self._variants:
            self._variants.append(molecule_variant)
            # Each variant might have its own conformers
            variant_key = molecule_variant.canonical_smiles(True)
            if variant_key is None:
                logger.warning(
                    f"Variant {molecule_variant.name} has no canonical SMILES, cannot manage conformers effectively."
                )
            else:
                self._conformer_sets[variant_key] = ConformerSet(
                    molecule_variant.rdkit_mol
                )
                logger.info(
                    f"Added variant {molecule_variant.name or molecule_variant.canonical_smiles(True)} to container."
                )

    def get_variants(self) -> list[Molecule]:
        """Returns all molecular variants in this container."""
        return self._variants

    def get_conformer_set(self, molecule_variant: Molecule) -> ConformerSet | None:
        """Retrieves the ConformerSet for a specific molecular variant."""
        return self._conformer_sets.get(molecule_variant.canonical_smiles(True) or "")

    def generate_all_conformers(
        self,
        num_conformers: int,
        rmsd_cutoff: float = 0.1,
        minimize: bool = True,
        use_random_coords: bool = False,
    ) -> None:
        """
        Generates, minimizes, and filters conformers for all variants in the container.
        """
        for variant in self._variants:
            conformer_set = self.get_conformer_set(variant)
            if conformer_set:
                logger.info(
                    f"Generating conformers for variant: {variant.name or variant.canonical_smiles(True)}"
                )
                conformer_set.generate_conformers(
                    num_conformers, use_random_coordinates=use_random_coords
                )
                if minimize:
                    conformer_set.minimize_all()
                conformer_set.eliminate_structurally_similar(rmsd_cutoff)
                logger.info(
                    f"Finished conformer generation for {variant.name or variant.canonical_smiles(True)}. Found {len(conformer_set.get_conformers())} unique conformers."
                )
            else:
                logger.warning(
                    f"No ConformerSet found for variant: {variant.name or variant.canonical_smiles(True)}. Skipping conformer generation for this variant."
                )

    def load_conformers_into_all_rdkit_mols(self) -> None:
        """Loads the generated conformers back into the RDKit Mol objects of all variants."""
        for variant in self._variants:
            conformer_set = self.get_conformer_set(variant)
            if conformer_set and conformer_set.get_conformers():
                conformer_set.load_into_rdkit_mol(variant.rdkit_mol)
            else:
                logger.info(
                    f"No conformers to load for variant: {variant.name or variant.canonical_smiles(True)}"
                )
