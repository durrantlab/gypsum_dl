"""Regression tests for the Durrant-lab filter substructure patterns.

The README is the only specification of this filter, so these tests pin the
documented list to the constant the filter actually compiles, and pin the
behavior of the two patterns that differ from the historical documentation.
"""

import os
import re

import pytest
from rdkit import Chem

from gypsum_dl.steps.smiles.DurrantLabFilter import (
    durrant_lab_contains_metal,
    metal_atomic_nums,
    metal_element_symbols,
    prohibited_smi_substrs_for_substruc,
)

README = os.path.join(os.path.dirname(os.path.dirname(__file__)), "README.md")

# Documented pattern that the code deliberately does not use: it requires an
# NH or NH2 nitrogen, so it misses N-substituted (internal) iminols.
TERMINAL_ONLY_IMINOL = "[$([NX2H1]),$([NX3H2])]=C[$([OH]),$([O-])]"


def _documented_patterns() -> list[str]:
    """Pull the filter patterns out of the README's Durrant-Lab Filters section.

    The list is scraped rather than duplicated here so that the test fails when
    the docs and the code drift apart, which is the failure mode it guards.

    Returns:
        The backtick-quoted pattern from each bullet in that section, in
        document order. Bullets without a backtick-quoted pattern (e.g.,
        "Metals") are skipped.
    """
    with open(README, encoding="utf-8") as f:
        text = f.read()

    section = re.split(r"^### Durrant-Lab Filters$", text, flags=re.MULTILINE)
    if len(section) < 2:
        pytest.fail("README.md has no '### Durrant-Lab Filters' section.")
    body = re.split(r"^#{2,3} ", section[1], flags=re.MULTILINE)[0]

    return re.findall(r"^- `([^`]+)`", body, flags=re.MULTILINE)


def test_readme_documents_the_patterns_the_code_uses() -> None:
    assert _documented_patterns() == prohibited_smi_substrs_for_substruc


@pytest.mark.parametrize("pattern", prohibited_smi_substrs_for_substruc)
def test_every_pattern_is_valid_smarts(pattern: str) -> None:
    assert Chem.MolFromSmarts(pattern) is not None


@pytest.mark.parametrize("pattern", prohibited_smi_substrs_for_substruc)
def test_no_pattern_depends_on_stereochemistry(pattern: str) -> None:
    # The filter runs before chiral and double-bond enumeration, which is only
    # sound while no pattern cares about stereochemistry. HasSubstructMatch
    # also ignores stereo unless asked, so a stereo pattern added here would
    # quietly match every stereoisomer rather than fail. Either way it needs
    # the step moved and useChirality=True, not just a new list entry.
    assert not any(marker in pattern for marker in "@/\\")


@pytest.mark.parametrize(
    ("pattern", "smiles", "expected"),
    [
        # Internal iminol: matched by the in-code pattern, missed by the
        # narrower pattern the README used to document.
        ("[$(N)]=C[$([OH]),$([O-])]", "CN=C(O)C", True),
        ("[$(N)]=C[$([OH]),$([O-])]", "N=C(O)C", True),
        ("[$(N)]=C[$([OH]),$([O-])]", "CN=C(OC)C", False),
        ("[$(N)]=C[$([OH]),$([O-])]", "CNC(=O)C", False),
        # Mistaken amide tautomer (enol-amide).
        ("[$(N)]C(=C)[$([OH]),$([O-])]", "NC(=C)O", True),
        ("[$(N)]C(=C)[$([OH]),$([O-])]", "CNC(=C)O", True),
        ("[$(N)]C(=C)[$([OH]),$([O-])]", "CC(N)=O", False),
    ],
)
def test_iminol_patterns_match_expected_molecules(
    pattern: str, smiles: str, expected: bool
) -> None:
    assert pattern in prohibited_smi_substrs_for_substruc
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    assert mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) is expected


def test_internal_iminol_is_only_caught_by_the_in_code_pattern() -> None:
    # Pins why the code pattern is broader than the previously documented one,
    # so a well-meaning "restore the README pattern" change fails here.
    mol = Chem.MolFromSmiles("CN=C(O)C")
    assert mol.HasSubstructMatch(Chem.MolFromSmarts(TERMINAL_ONLY_IMINOL)) is False
    assert TERMINAL_ONLY_IMINOL not in prohibited_smi_substrs_for_substruc


# The eleven metals the filter rejected before the list was completed.
ORIGINAL_METALS = ["Al", "V", "Fe", "Co", "Cu", "Zn", "Mo", "Cd", "Au", "Pb", "Bi"]

# Metals that used to pass the filter, which is the gap these tests close.
PREVIOUSLY_PERMITTED_METALS = [
    "Mg",
    "Mn",
    "Ni",
    "Hg",
    "Pt",
    "Ag",
    "Pd",
    "Ru",
    "Sn",
    "Ti",
    "Cr",
]


@pytest.mark.parametrize("symbol", metal_element_symbols)
def test_every_metal_symbol_is_a_real_element(symbol: str) -> None:
    # A typo in the symbol list would silently stop filtering that metal, since
    # a substring that matches nothing is indistinguishable from one that has
    # nothing to match.
    mol = Chem.MolFromSmiles(f"[{symbol}]")
    assert mol is not None
    assert mol.GetAtomWithIdx(0).GetSymbol() == symbol
    assert mol.GetAtomWithIdx(0).GetAtomicNum() <= 92


def _contains_metal(smiles: str) -> bool:
    """Run the metal check on a SMILES string.

    The check takes a molecule, so each test would otherwise repeat the parse
    and the guard against a SMILES that RDKit rejects.

    Args:
        smiles: The SMILES string to test.

    Returns:
        Whether the parsed molecule contains a metal atom.
    """
    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None
    return durrant_lab_contains_metal(mol)


@pytest.mark.parametrize("symbol", metal_element_symbols)
def test_every_metal_symbol_is_filtered_by_atomic_number(symbol: str) -> None:
    atom = Chem.MolFromSmiles(f"[{symbol}]").GetAtomWithIdx(0)
    assert atom.GetAtomicNum() in metal_atomic_nums


@pytest.mark.parametrize("symbol", ORIGINAL_METALS + PREVIOUSLY_PERMITTED_METALS)
def test_metal_containing_smiles_are_rejected(symbol: str) -> None:
    assert _contains_metal(f"[{symbol}]") is True
    assert _contains_metal(f"CC(=O)[O-].[{symbol}+2]") is True


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=O)Oc1ccccc1C(=O)O",  # Aspirin.
        "C[C@H](N)C(=O)O",  # Bracketed chiral carbon, not cadmium or cobalt.
        "c1cc[nH]c1",  # Bracketed aromatic nitrogen, not sodium.
        "C[N+](C)(C)C",  # Charged nitrogen, not niobium.
        "[I-]",  # Iodide, not indium.
        "O=[PH](O)O",  # Bracketed phosphorus, not lead or platinum.
        "[Si](C)(C)C",  # A metalloid, left alone on purpose.
        "[2H]C(Cl)(Cl)Cl",  # Isotope label on a non-metal.
    ],
)
def test_non_metal_bracket_atoms_are_kept(smiles: str) -> None:
    # Bracket atoms whose symbols share a first letter with a metal. The old
    # SMILES substring test had to be written carefully to leave these alone.
    assert _contains_metal(smiles) is False


def test_krypton_is_not_mistaken_for_potassium() -> None:
    # Regression: the substring "[K" matched krypton as well as potassium.
    assert _contains_metal("[Kr]") is False
    assert _contains_metal("[K+]") is True


@pytest.mark.parametrize("smiles", ["[23Na+]", "CC(=O)[O-].[64Cu+2]"])
def test_isotope_labeled_metals_are_rejected(smiles: str) -> None:
    # Regression: an isotope label puts digits between the bracket and the
    # symbol, so the substring "[Na" never matched "[23Na+]".
    assert _contains_metal(smiles) is True
