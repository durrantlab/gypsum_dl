"""
This module contains classes and functions for processing individual molecules
(variants).
"""

from typing import Any

import copy
from collections.abc import Container

from loguru import logger
from molvs import standardize_smiles as ssmiles
from rdkit import Chem, RDLogger
from rdkit.Chem import BondStereo

from gypsum_dl import handlers

RDLogger.DisableLog("rdApp.*")  # type: ignore


class Molecule:
    """
    A class that wraps around an RDKit Mol object, focusing on its core
    structural and chemical properties.
    """

    def __init__(self, rdkit_mol: Chem.Mol, name: str = ""):
        """
        Initialize the Molecule object with an RDKit Mol object.
        Sanitization is performed upon initialization.
        """
        if rdkit_mol is None:
            raise ValueError("RDKit Mol object cannot be None.")

        # Ensure the RDKit Mol is properly sanitized
        self._rdkit_mol = handlers.check_sanitization(rdkit_mol)
        if self._rdkit_mol is None:
            raise ValueError("Could not sanitize the provided RDKit Mol.")

        self._name = name
        self._orig_smiles = Chem.MolToSmiles(
            self._rdkit_mol, isomericSmiles=True, canonical=False
        )
        self._can_smiles: str | None = None
        self._can_smiles_noh: str | None = None
        self._standardized_smiles: str | None = None
        self._mol_props: dict[str, Any] = {}
        self._genealogy: list[str] = []

        # Cached properties
        self._nonaro_ring_atom_idx: list[list[int]] | None = None
        self._chiral_centers_assigned: list[Any] | None = None
        self._chiral_centers_unassigned: list[Any] | None = None
        self._has_bizarre_substructure: bool | None = None
        self._fragments: Container[Chem.Mol] | None = None

    @classmethod
    def from_smiles(cls, smiles: str, name: str = "") -> "Molecule":
        """Factory method to create a Molecule from a SMILES string."""
        try:
            # sanitize=False to respect double-bond stereochemistry initially
            mol = Chem.MolFromSmiles(smiles, sanitize=False)
            if mol is None:
                raise ValueError(f"Could not create RDKit Mol from SMILES: {smiles}")
            return cls(mol, name)
        except Exception as e:
            logger.error(f"Failed to create Molecule from SMILES '{smiles}': {e}")
            raise

    @classmethod
    def from_rdkit_mol(cls, rdkit_mol: Chem.Mol, name: str = "") -> "Molecule":
        """Factory method to create a Molecule from an existing RDKit Mol object."""
        return cls(rdkit_mol, name)

    @property
    def rdkit_mol(self) -> Chem.Mol:
        return self._rdkit_mol

    @property
    def name(self) -> str:
        return self._name

    @name.setter
    def name(self, new_name: str) -> None:
        self._name = new_name

    @property
    def original_smiles(self) -> str:
        return self._orig_smiles

    def canonical_smiles(self, include_hydrogens: bool = True) -> str | None:
        """
        Get the canonical SMILES string associated with this object.
        Caches the result.
        """
        if include_hydrogens:
            if self._can_smiles is None:
                try:
                    self._can_smiles = Chem.MolToSmiles(
                        self._rdkit_mol, isomericSmiles=True, canonical=True
                    )
                except Exception:
                    logger.warning(
                        f"Couldn't generate canonical SMILES for {self._orig_smiles} ({self._name})."
                    )
                    self._can_smiles = None
            return self._can_smiles
        else:
            if self._can_smiles_noh is None:
                # Temporarily deprotonate to get SMILES without explicit hydrogens
                deprotonated_mol = handlers.try_deprotanation(
                    copy.deepcopy(self._rdkit_mol)
                )
                if deprotonated_mol is not None:
                    self._can_smiles_noh = Chem.MolToSmiles(
                        deprotonated_mol, isomericSmiles=True, canonical=True
                    )
                else:
                    logger.warning(
                        f"Could not deprotonate molecule for non-H canonical SMILES: {self._orig_smiles}"
                    )
                    self._can_smiles_noh = None
            return self._can_smiles_noh

    @property
    def standardized_smiles(self) -> str | None:
        """Returns the standardized SMILES string, computing and caching if necessary."""
        if self._standardized_smiles is None:
            # It's important to pass a string to ssmiles.
            can_smi = self.canonical_smiles(include_hydrogens=True)
            if can_smi is None:
                logger.warning(
                    f"Cannot standardize SMILES for {self.original_smiles} because canonical SMILES is None."
                )
                self._standardized_smiles = None
            else:
                try:
                    self._standardized_smiles = ssmiles(can_smi)  # type: ignore
                except Exception:
                    logger.info(
                        f"Could not standardize SMILES for {can_smi}. Skipping."
                    )
                    self._standardized_smiles = (
                        can_smi  # Fallback to canonical if standardization fails
                    )
        return self._standardized_smiles

    def add_to_genealogy(self, step_description: str) -> None:
        """Adds a step to the molecule's genealogy (preparation history)."""
        self._genealogy.append(step_description)

    @property
    def genealogy(self) -> list[str]:
        return self._genealogy

    def set_property(self, key: str, value: Any) -> None:
        """Set a custom molecular property."""
        self._mol_props[key] = value
        self.update_rdkit_mol_properties()  # Update RDKit Mol too

    def get_property(self, key: str) -> Any:
        """Get a custom molecular property."""
        return self._mol_props.get(key)

    def update_rdkit_mol_properties(self) -> None:
        """Copies all stored properties to the RDKit Mol object."""
        self._rdkit_mol.SetProp("SMILES", self.canonical_smiles(True) or "N/A")
        if self.original_smiles:
            self._rdkit_mol.SetProp("ORIGINAL_SMILES", self.original_smiles)
        if self._name:
            self._rdkit_mol.SetProp("_Name", self._name)
        if self.genealogy:
            self._rdkit_mol.SetProp("Genealogy", "\n".join(self.genealogy))
        for prop_key, prop_val in self._mol_props.items():
            try:
                self._rdkit_mol.SetProp(prop_key, str(prop_val))
            except Exception as e:
                logger.warning(
                    f"Could not set RDKit property '{prop_key}' with value '{prop_val}': {e}"
                )

    def add_explicit_hydrogens(self) -> None:
        """Adds explicit hydrogen atoms to the RDKit Mol object."""
        new_mol = handlers.try_reprotanation(self._rdkit_mol)
        if new_mol is None:
            raise RuntimeError("Failed to add explicit hydrogens to molecule.")
        self._rdkit_mol = new_mol
        self._can_smiles = None  # Invalidate cached SMILES
        self._can_smiles_noh = None

    def remove_explicit_hydrogens(self) -> None:
        """Removes explicit hydrogen atoms from the RDKit Mol object."""
        new_mol = handlers.try_deprotanation(self._rdkit_mol)
        if new_mol is None:
            raise RuntimeError("Failed to remove explicit hydrogens from molecule.")
        self._rdkit_mol = new_mol
        self._can_smiles = None  # Invalidate cached SMILES
        self._can_smiles_noh = None

    @property
    def non_aromatic_ring_atom_indices(self) -> list[list[int]]:
        """Identifies and returns indices of atoms in non-aromatic rings, caching the result."""
        if self._nonaro_ring_atom_idx is None:
            if self._rdkit_mol is None:
                self._nonaro_ring_atom_idx = []
            else:
                ssr = Chem.GetSymmSSSR(self._rdkit_mol)
                ring_indices = [list(ssr[i]) for i in range(len(ssr))]
                non_aromatic_rings = []
                for ring_idx_set in ring_indices:
                    if any(
                        not self._rdkit_mol.GetAtomWithIdx(atm_idx).GetIsAromatic()
                        for atm_idx in ring_idx_set
                    ):
                        non_aromatic_rings.append(ring_idx_set)
                self._nonaro_ring_atom_idx = non_aromatic_rings
        return self._nonaro_ring_atom_idx

    @property
    def chiral_centers_assigned(self) -> list[Any]:
        """Returns a list of assigned chiral centers, caching the result."""
        if self._chiral_centers_assigned is None:
            if self._rdkit_mol is None:
                self._chiral_centers_assigned = []
            else:
                self._chiral_centers_assigned = Chem.FindMolChiralCenters(
                    self._rdkit_mol, includeUnassigned=False
                )
        return self._chiral_centers_assigned

    @property
    def chiral_centers_unassigned(self) -> list[Any]:
        """Returns a list of chiral centers including unassigned ones, caching the result."""
        if self._chiral_centers_unassigned is None:
            if self._rdkit_mol is None:
                self._chiral_centers_unassigned = []
            else:
                self._chiral_centers_unassigned = Chem.FindMolChiralCenters(
                    self._rdkit_mol, includeUnassigned=True
                )
        return self._chiral_centers_unassigned

    @property
    def double_bonds_without_stereochemistry(self) -> list[int]:
        """Returns indices of double bonds without specified stereochemistry."""
        if self._rdkit_mol is None:
            return []
        return [
            b.GetIdx()
            for b in self._rdkit_mol.GetBonds()
            if b.GetBondTypeAsDouble() == 2 and b.GetStereo() is BondStereo.STEREONONE
        ]

    @property
    def has_bizarre_substructure(self) -> bool:
        """
        Checks for and caches whether the molecule contains improbable substructures.
        """
        if self._has_bizarre_substructure is None:
            if self._rdkit_mol is None:
                self._has_bizarre_substructure = True  # No mol is bizarre
                return True

            prohibited_substructures = [
                "O(=*)-*",
                "C(=[CH2])[OH]",
                "C(=[CH2])[O-]",
                "C=C([OH])[OH]",
                "C=C([O-])[OH]",
                "C=C([O-])[O-]",
                "[C-]",
                "[c-]",
            ]

            for s in prohibited_substructures:
                pattrn = Chem.MolFromSmarts(s)
                if pattrn is not None and self._rdkit_mol.HasSubstructMatch(pattrn):
                    logger.info(
                        f"Detected unusual substructure: {s} in {self.canonical_smiles(True)}"
                    )
                    self._has_bizarre_substructure = True
                    return True
            self._has_bizarre_substructure = False
        return self._has_bizarre_substructure

    @property
    def fragments(self) -> Container[Chem.Mol]:
        """Divides the molecule into fragments if disconnected, caching the result."""
        if self._fragments is None:
            if "." not in self._orig_smiles:
                self._fragments = [self.rdkit_mol]  # A list containing only itself
            else:
                self._fragments = Chem.GetMolFrags(self._rdkit_mol, asMols=True)
        return self._fragments

    def count_hydrogens_bound_to_carbons(self) -> int:
        """Count the number of Hydrogens bound to carbons."""
        if self._rdkit_mol is None:
            return 0

        total_hydrogens_counted = 0
        for atom in self._rdkit_mol.GetAtoms():
            if atom.GetSymbol() == "C":
                total_hydrogens_counted += atom.GetTotalNumHs(includeNeighbors=True)
        return total_hydrogens_counted

    def __hash__(self) -> int:
        """Allows hashing based on the canonical SMILES."""
        can_smi = self.canonical_smiles(include_hydrogens=True)
        return hash(can_smi) if can_smi else 0

    def __eq__(self, other: object) -> bool:
        """Compares Molecule objects based on canonical SMILES."""
        if not isinstance(other, Molecule):
            return NotImplemented
        return self.canonical_smiles(True) == other.canonical_smiles(True)

    def __ne__(self, other: object) -> bool:
        return not self.__eq__(other)

    def __lt__(self, other: "Molecule") -> bool:
        if not isinstance(other, Molecule):
            return NotImplemented
        return (self.canonical_smiles(True) or "") < (
            other.canonical_smiles(True) or ""
        )

    def __le__(self, other: "Molecule") -> bool:
        if not isinstance(other, Molecule):
            return NotImplemented
        return (self.canonical_smiles(True) or "") <= (
            other.canonical_smiles(True) or ""
        )

    def __gt__(self, other: "Molecule") -> bool:
        if not isinstance(other, Molecule):
            return NotImplemented
        return (self.canonical_smiles(True) or "") > (
            other.canonical_smiles(True) or ""
        )

    def __ge__(self, other: "Molecule") -> bool:
        if not isinstance(other, Molecule):
            return NotImplemented
        return (self.canonical_smiles(True) or "") >= (
            other.canonical_smiles(True) or ""
        )
