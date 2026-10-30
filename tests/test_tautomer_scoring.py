"""Tests for ranking tautomers by RDKit's tautomer score.

Through 2.0.0 the tautomer step kept the tautomers with the lowest UFF energy
on an unminimized embedding, which cannot rank tautomers, so the kept form of
an enolizable compound was close to arbitrary and changed with the seed.
"""

import pytest

from gypsum_dl import MyMol, chem_utils
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.parallelizer import seed_generators
from gypsum_dl.steps.smiles.MakeTautomers import make_tauts


def _mols(*smiles: str) -> list[MyMol.MyMol]:
    """Build one candidate per SMILES.

    Args:
        *smiles: The candidates, in the order the selector should see them.

    Returns:
        The candidates as MyMol objects.
    """
    return [MyMol.MyMol(s) for s in smiles]


def _canon(smiles: str) -> str:
    """Express a SMILES the way the selector reports it.

    Args:
        smiles: Any valid SMILES.

    Returns:
        The SMILES MyMol produces for it.
    """
    return MyMol.MyMol(smiles).smiles()


@pytest.mark.parametrize(
    ("favored", "other"),
    [
        pytest.param("CC(C)=O", "C=C(C)O", id="ketone_over_enol"),
        pytest.param("CCOC(C)=O", "C=C(O)OCC", id="ester_over_enol"),
        pytest.param(
            "O=c1cc[nH]c2ccccc12", "Oc1ccnc2ccccc12", id="quinolone_over_hydroxy"
        ),
    ],
)
def test_scorer_prefers_the_conventional_tautomer(favored: str, other: str) -> None:
    # The disfavored form is listed first and named as the parent, so neither
    # input order nor the tie-break can be what picks the favored one.
    kept = chem_utils.pick_highest_scoring_tautomers(_mols(other, favored), 1, {other})

    assert [m.smiles() for m in kept] == [_canon(favored)]


@pytest.mark.parametrize(
    "parent",
    ["CCOC(=O)Nc1nnn[nH]1", "CCOC(=O)Nc1nn[nH]n1"],
    ids=["1H", "2H"],
)
def test_ties_go_to_the_form_the_step_started_from(parent: str) -> None:
    one_h, two_h = _mols("CCOC(=O)Nc1nnn[nH]1", "CCOC(=O)Nc1nn[nH]n1")
    assert chem_utils.tautomer_score(one_h) == chem_utils.tautomer_score(two_h)

    kept = chem_utils.pick_highest_scoring_tautomers([one_h, two_h], 1, {parent})

    assert [m.smiles() for m in kept] == [_canon(parent)]


def test_each_protonation_state_gets_a_slot_before_any_gets_two() -> None:
    # Scores do not compare across protonation states, so the selector keeps
    # the round-robin across formula-and-charge groups that the energy ranking
    # uses: with two slots, the cation takes the second rather than the enol.
    kept = chem_utils.pick_highest_scoring_tautomers(
        _mols("C=C(C)O", "CC(C)=O", "C[NH3+]"), 2, set()
    )

    assert [m.smiles() for m in kept] == [_canon("CC(C)=O"), _canon("C[NH3+]")]


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
@pytest.mark.parametrize(
    "smiles", ["CCOC(C)=O", "CCC(=O)CC(C)(C)C", "CC(C)CC(=O)OCC(C)(C)C"]
)
def test_make_tauts_keeps_the_conventional_form_whatever_the_seed(
    smiles: str, seed: int
) -> None:
    # The last two are compounds 896 and 937 of RDKit's bundled NCI sample,
    # which 2.0.0 kept as enols under some seeds and as the keto or ester
    # form under others.
    seed_generators(seed)
    contnr = MolContainer(smiles, "compound", 0, {})
    contnr.add_smiles(smiles)

    make_tauts([contnr], 1, 3, 1, "serial", False, None)

    assert [m.smiles() for m in contnr.mols] == [_canon(smiles)]
    # Scoring needs no throwaway conformers, so none were embedded.
    assert contnr.probe_energies == {}
