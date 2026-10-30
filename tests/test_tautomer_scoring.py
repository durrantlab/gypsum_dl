"""Tests for ranking tautomers by RDKit's tautomer score.

Through 2.0.0 the tautomer step kept the tautomers with the lowest UFF energy
on an unminimized embedding, which cannot rank tautomers, so the kept form of
an enolizable compound was close to arbitrary and changed with the seed.
"""

import pytest

from molvs import tautomer as molvs_tautomer

from gypsum_dl import MyMol, chem_utils
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.parallelizer import seed_generators
from gypsum_dl.steps.smiles.MakeTautomers import (
    IMINE_FREE_TRANSFORMS,
    make_tauts,
    parallel_make_taut,
)

ENAMINO_ESTER = "CCOC(=O)C=C(C)N"
ENAMINO_ESTER_IMINE = "CCOC(=O)CC(C)=N"


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


def test_scorer_alone_still_prefers_hydroximic_acids() -> None:
    # Documents why _HYDROXIMIC_ACID exists. If a future RDKit stops awarding
    # the oxime bonus that causes this, the rule may no longer be needed.
    hydroxamic, hydroximic = _mols("CC(=O)NO", "CC(O)=NO")

    assert chem_utils.tautomer_score(hydroximic) > chem_utils.tautomer_score(hydroxamic)


@pytest.mark.parametrize(
    ("hydroxamic", "hydroximic"),
    [
        pytest.param("CC(=O)NO", "CC(O)=NO", id="neutral"),
        pytest.param("CC(=O)N[O-]", "CC([O-])=NO", id="anion"),
    ],
)
def test_hydroxamic_acids_rank_above_their_hydroximic_tautomers(
    hydroxamic: str, hydroximic: str
) -> None:
    # The hydroximic form is listed first and named as the parent, so only the
    # rule can be what keeps the hydroxamic acid.
    kept = chem_utils.pick_highest_scoring_tautomers(
        _mols(hydroximic, hydroxamic), 1, {hydroximic}
    )

    assert [m.smiles() for m in kept] == [_canon(hydroxamic)]


@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param("NC(Cc1ccc(O)cc1)C(=O)NO", id="nci_667"),
        pytest.param("O=C(CCCCCCC(=O)Nc1ccccc1)NO", id="vorinostat"),
    ],
)
def test_make_tauts_keeps_hydroxamic_acids_without_the_filters(smiles: str) -> None:
    # With --use_durrant_lab_filters on, the iminol pattern already rejects the
    # hydroximic form; this covers runs with the filters off.
    contnr = MolContainer(smiles, "compound", 0, {})
    contnr.add_smiles(smiles)

    make_tauts([contnr], 1, 3, 1, "serial", False, None)

    assert [m.smiles() for m in contnr.mols] == [_canon(smiles)]


def test_nitronic_acids_rank_below_their_nitro_form() -> None:
    # The nitronic form is listed first and named as the parent, so only the
    # rule can be what keeps the nitro group.
    kept = chem_utils.pick_highest_scoring_tautomers(
        _mols("CC=[N+]([O-])O", "CC[N+](=O)[O-]"), 1, {"CC=[N+]([O-])O"}
    )

    assert [m.smiles() for m in kept] == [_canon("CC[N+](=O)[O-]")]


def test_imine_free_transforms_drop_exactly_the_imine_shifts() -> None:
    # The transforms are selected by name, so a MolVS release that renamed or
    # added one would otherwise change the rule silently.
    kept = {t.name for t in IMINE_FREE_TRANSFORMS}
    dropped = {t.name for t in molvs_tautomer.TAUTOMER_TRANSFORMS} - kept

    assert dropped == {
        "aliphatic imine f",
        "aliphatic imine r",
        "special imine f",
        "special imine r",
    }


def test_worker_flags_only_the_forms_that_need_an_imine_shift() -> None:
    contnr = MolContainer(ENAMINO_ESTER, "enamino_ester", 0, {})
    contnr.add_smiles(ENAMINO_ESTER)

    tauts = parallel_make_taut(contnr.mols[0], contnr.contnr_props(), 3, True)

    flags = {
        chem_utils.constitution_smiles(t.smiles()): t.needs_imine_shift
        for t in tauts or []
    }
    assert flags[chem_utils.constitution_smiles(ENAMINO_ESTER)] is False
    assert flags[chem_utils.constitution_smiles(ENAMINO_ESTER_IMINE)] is True


def test_imine_shift_forms_rank_below_the_rest_but_keep_a_later_slot() -> None:
    # The imine is named as the parent, so only the flag can put the enamine
    # first; with two slots, the imine is still offered as the second variant.
    enamine, imine = _mols(ENAMINO_ESTER, ENAMINO_ESTER_IMINE)
    imine.needs_imine_shift = True

    one = chem_utils.pick_highest_scoring_tautomers(
        [imine, enamine], 1, {ENAMINO_ESTER_IMINE}
    )
    two = chem_utils.pick_highest_scoring_tautomers(
        [imine, enamine], 2, {ENAMINO_ESTER_IMINE}
    )

    assert [m.smiles() for m in one] == [_canon(ENAMINO_ESTER)]
    assert [m.smiles() for m in two] == [
        _canon(ENAMINO_ESTER),
        _canon(ENAMINO_ESTER_IMINE),
    ]


def _kept_constitution(smiles: str) -> str | None:
    """Run the tautomer step at production settings and report what it kept.

    Args:
        smiles: The input form.

    Returns:
        The constitution of the one tautomer kept.
    """
    contnr = MolContainer(smiles, "compound", 0, {})
    contnr.add_smiles(smiles)
    make_tauts([contnr], 1, 3, 1, "serial", False, None)
    assert len(contnr.mols) == 1
    return chem_utils.constitution_smiles(contnr.mols[0].smiles())


# Compounds from RDKit's bundled NCI sample (and nizatidine, CHEMBL3183075)
# whose conjugated enamine the score turned into an imine.
@pytest.mark.parametrize(
    "smiles",
    [
        pytest.param(ENAMINO_ESTER, id="nci_1086_enamino_ester"),
        pytest.param("COc1ccc(C=C2NC(=O)NC2=O)cc1", id="nci_1232_arylidene"),
        pytest.param("O=C1NC(=S)NC1=Cc1ccccc1O", id="nci_2592_arylidene"),
        pytest.param("CCC(N)=CC(=O)c1ccccc1", id="nci_3153_enaminone"),
        pytest.param("CC1=CC(C)(C)NC(=S)N1", id="nci_3270_cyclic"),
        pytest.param("CNC(=C[N+](=O)[O-])NCCSCc1csc(CN(C)C)n1", id="nizatidine"),
    ],
)
def test_make_tauts_keeps_conjugated_enamines(smiles: str) -> None:
    assert _kept_constitution(smiles) == chem_utils.constitution_smiles(smiles)


# Changes the same sample showed the score gets right, which the imine and
# nitronic acid rules must leave alone.
@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        pytest.param(
            "N=c1[nH]cc([N+](=O)[O-])s1", "Nc1ncc([N+](=O)[O-])s1", id="aminothiazole"
        ),
        pytest.param(
            "CC(=O)Nc1ccc2c(O)ccnc2n1",
            "CC(=O)Nc1ccc2c(=O)cc[nH]c2n1",
            id="lactam",
        ),
        pytest.param("Nc1nc(S)c2nc[nH]c2n1", "Nc1nc(=S)c2nc[nH]c2[nH]1", id="thione"),
    ],
)
def test_make_tauts_still_corrects_minor_input_forms(
    smiles: str, expected: str
) -> None:
    assert _kept_constitution(smiles) == chem_utils.constitution_smiles(expected)
