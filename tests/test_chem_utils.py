"""Unit tests for the variant-selection helpers in chem_utils."""

import pytest

from gypsum_dl import MyMol, chem_utils
from gypsum_dl.MolContainer import MolContainer


def test_uniq_mols_in_list_removes_duplicates() -> None:
    mols = [MyMol.MyMol("CCO"), MyMol.MyMol("OCC"), MyMol.MyMol("CCC")]
    assert len(chem_utils.uniq_mols_in_list(mols)) == 2


def test_uniq_mols_in_list_keeps_molecules_with_unknown_smiles() -> None:
    # Regression: the seen-set was keyed on smiles() directly, and smiles()
    # reports failure as None. Two distinct molecules that both failed to
    # canonicalize therefore shared the key None, so the second was discarded
    # as a duplicate of the first.
    first = MyMol.MyMol("CCO")
    second = MyMol.MyMol("CCCCCC")
    first.can_smi = None
    second.can_smi = None

    kept = chem_utils.uniq_mols_in_list([first, second])

    assert len(kept) == 2


def test_remove_highly_charged_molecules_tolerates_unknown_smiles() -> None:
    # The discard warning concatenated smiles() into a string, which is a
    # TypeError in the main process once canonicalization has failed.
    neutral = MyMol.MyMol("CCO")
    charged = MyMol.MyMol("[NH4+].[NH4+].[NH4+].[NH4+].[NH4+]")
    charged.can_smi = None

    kept = chem_utils.remove_highly_charged_molecules([neutral, charged])

    assert len(kept) == 1


def test_remove_highly_charged_molecules_discards_outliers() -> None:
    neutral = MyMol.MyMol("CCO")
    charged = MyMol.MyMol("[NH4+].[NH4+].[NH4+].[NH4+].[NH4+]")
    kept = chem_utils.remove_highly_charged_molecules([neutral, charged])
    assert len(kept) == 1
    assert kept[0].smiles() == neutral.smiles()


def test_remove_highly_charged_molecules_keeps_similar_charges() -> None:
    mols = [MyMol.MyMol("CCO"), MyMol.MyMol("[NH4+]")]
    assert len(chem_utils.remove_highly_charged_molecules(mols)) == 2


def test_pick_lowest_enrgy_mols_returns_all_when_under_limit() -> None:
    mols = [MyMol.MyMol("CCO"), MyMol.MyMol("CCC")]
    assert len(chem_utils.pick_lowest_enrgy_mols(mols, 5, 1)) == 2


def test_pick_lowest_enrgy_mols_returns_lowest_energy_mol() -> None:
    """Regression: indices must come from mols_3d, not mol_lst (which has a
    different order after random_sample / set dedup)."""

    class _FakeConf:
        def __init__(self, energy):
            self.energy = energy

    smiles = ["C", "CC", "CCC"]
    mols = [MyMol.MyMol(s) for s in smiles]

    # Assign stable energies and pre-populate conformers so no 3D embedding runs.
    energies = {"C": 100.0, "CC": 50.0, "CCC": 10.0}
    for mol in mols:
        e = energies[mol.smiles()]
        mol.conformers = [_FakeConf(e)]
        mol.make_first_3d_conf_no_min = (
            lambda m=mol: None
        )  # no-op; conformers already set

    # thoroughness=3 so random_sample draws all 3 candidates; run 20 times to
    # defeat the internal shuffle.
    for _ in range(20):
        kept = chem_utils.pick_lowest_enrgy_mols(mols, 1, 3)
        assert len(kept) == 1
        assert (
            kept[0].smiles() == "CCC"
        ), f"Expected lowest-energy mol 'CCC' but got '{kept[0].smiles()}'"


def test_pick_lowest_enrgy_mols_leaves_candidates_unchanged() -> None:
    # Ranking used to build the ranking conformer on the candidates
    # themselves, which reprotonates rdkit_mol and attaches 3D coordinates. A
    # variant that survived a pruning step then differed from a variant whose
    # container never needed pruning, which surfaced downstream as real 3D
    # coordinates in the SDF of a 2d_output_only run.
    mols = [MyMol.MyMol(s) for s in ["CCO", "CCC", "CCCC"]]
    atom_counts = [m.rdkit_mol.GetNumAtoms() for m in mols]

    # thoroughness=3 so every candidate gets ranked, not just a sample of one.
    kept = chem_utils.pick_lowest_enrgy_mols(mols, 1, 3)

    assert len(kept) == 1
    assert any(kept[0] is m for m in mols)

    for mol, atom_count in zip(mols, atom_counts):
        assert mol.conformers == []
        assert mol.rdkit_mol.GetNumAtoms() == atom_count
        assert mol.rdkit_mol.GetNumConformers() == 0


def test_bst_for_each_contnr_no_opt_repopulates_containers() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    mol = MyMol.MyMol("CCO")
    mol.contnr_idx = 0
    chem_utils.bst_for_each_contnr_no_opt([contnr], [mol], 1, 1)
    assert len(contnr.mols) == 1


def test_bst_for_each_contnr_no_opt_carries_over_originals() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    chem_utils.bst_for_each_contnr_no_opt([contnr], [], 1, 1)
    assert len(contnr.mols) == 1


def test_bst_for_each_contnr_no_opt_can_discard_originals() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    chem_utils.bst_for_each_contnr_no_opt(
        [contnr], [], 1, 1, crry_ovr_frm_lst_step_if_no_fnd=False
    )
    assert contnr.mols == []


def test_bst_for_each_contnr_no_opt_rejects_a_foreign_container_index() -> None:
    # A candidate tagged for a container missing from the list used to be
    # dropped without a word, leaving its real container on the previous
    # step's variants.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    mol = MyMol.MyMol("CCO")
    mol.contnr_idx = 1

    with pytest.raises(Exception, match="contnr_idx"):
        chem_utils.bst_for_each_contnr_no_opt([contnr], [mol], 1, 1)


def test_bst_for_each_contnr_no_opt_names_the_step_in_its_warning(capsys) -> None:
    # The warning spoke of low-energy conformations whichever step called it,
    # although every SMILES-stage caller works before any conformer exists.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")

    chem_utils.bst_for_each_contnr_no_opt([contnr], [], 1, 1, variant_desc="tautomers")

    log = " ".join(capsys.readouterr().out.split())
    assert "No tautomers remained for CCO (ethanol)" in log
    assert "conformations" not in log


def test_pick_lowest_enrgy_mols_dedups_in_first_seen_order() -> None:
    # Regression: deduplication used list(set(...)), whose order follows
    # PYTHONHASHSEED because MyMol hashes its canonical SMILES string. The
    # deduplicated list feeds random_sample, so a fixed random_seed did not
    # pin down which variants advanced.
    mols = [MyMol.MyMol(s) for s in ["CCO", "CCC", "CCCC", "CCCCC", "OCC"]]

    # Limit above the number of distinct molecules, so the deduplicated list is
    # returned as is and its order is what gets asserted.
    kept = chem_utils.pick_lowest_enrgy_mols(mols, 10, 1)

    assert [m.smiles() for m in kept] == [m.smiles() for m in mols[:4]]


def test_remove_highly_charged_molecules_is_order_independent() -> None:
    # Regression: the reference charge came from abs_charges.index(min(...)),
    # which breaks ties by list position. With a -1 and a +1 form both one unit
    # from neutral, the signed reference charge (and so which forms survived)
    # followed whatever order the variants happened to arrive in.
    minus_one = MyMol.MyMol("CC(=O)[O-]")
    plus_one = MyMol.MyMol("CC[NH3+]")
    plus_five = MyMol.MyMol("[NH4+].[NH4+].[NH4+].[NH4+].[NH4+]")

    forward = chem_utils.remove_highly_charged_molecules(
        [plus_one, minus_one, plus_five]
    )
    reverse = chem_utils.remove_highly_charged_molecules(
        [minus_one, plus_one, plus_five]
    )

    assert {m.smiles() for m in forward} == {m.smiles() for m in reverse}
    # The reference is the -1 form either way, which puts the +5 form six units
    # away and so outside the window.
    assert len(forward) == 2


def test_pick_lowest_enrgy_mols_ranks_a_failed_energy_last() -> None:
    # Regression: a conformer whose force field failed carried the sentinel
    # 9999, which is not an upper bound on UFF energy. A strained or large
    # ligand above that value lost to a variant that was never scored at all,
    # and the SDF then reported 9999 as if it were a measurement.
    class _FakeConf:
        def __init__(self, energy: float) -> None:
            self.energy = energy

    strained = MyMol.MyMol("CCO")
    strained.conformers = [_FakeConf(12000.0)]

    unscored = MyMol.MyMol("CCC")
    unscored.conformers = [_FakeConf(float("inf"))]

    # thoroughness=2 so both candidates are ranked; repeat to defeat the
    # shuffle inside random_sample.
    for _ in range(20):
        kept = chem_utils.pick_lowest_enrgy_mols([strained, unscored], 1, 2)
        assert len(kept) == 1
        assert kept[0].smiles() == strained.smiles()


def test_first_conf_energy_reuses_a_cached_probe_energy() -> None:
    # The ranking embedding runs in the dispatching process, not through the
    # parallelizer, and a variant that survives one pruning step is a candidate
    # again at the next one. Without the cache the same structure is embedded
    # once per SMILES step.
    mol = MyMol.MyMol("CCO")
    cache: dict[str, float | None] = {}

    first = chem_utils.first_conf_energy(mol, cache)

    assert cache == {mol.smiles(): first}

    embedded = False

    def refuse() -> None:
        """Fail the test if a second embedding is attempted.

        Raises:
            AssertionError: Always, since a cache hit must not embed.
        """
        nonlocal embedded
        embedded = True
        raise AssertionError("embedded again despite a cache hit")

    probe = MyMol.MyMol("CCO")
    probe.make_first_3d_conf_no_min = refuse

    assert chem_utils.first_conf_energy(probe, cache) == first
    assert not embedded


def test_first_conf_energy_does_not_cache_an_existing_conformer_energy() -> None:
    # An energy read off coordinates the molecule already carries describes
    # that coordinate set, not the structure, so it must not stand in for the
    # ranking energy of another copy of the same structure.
    class _FakeConf:
        def __init__(self, energy: float) -> None:
            self.energy = energy

    mol = MyMol.MyMol("CCO")
    mol.conformers = [_FakeConf(42.0)]
    cache: dict[str, float | None] = {}

    assert chem_utils.first_conf_energy(mol, cache) == 42.0
    assert cache == {}


def test_first_conf_energy_skips_the_cache_for_an_unknown_smiles() -> None:
    # smiles() reports a failed canonicalization as None, which is no kind of
    # cache key: two unrelated molecules would share it.
    mol = MyMol.MyMol("CCO")
    mol.can_smi = None
    cache: dict[str, float | None] = {}

    chem_utils.first_conf_energy(mol, cache)

    assert cache == {}


def test_bst_for_each_contnr_no_opt_records_ranking_energies_on_the_container() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    mols = []
    for smi in ["CCO", "CCCO", "CCCCO"]:
        mol = MyMol.MyMol(smi)
        mol.contnr_idx = 0
        mols.append(mol)

    assert contnr.probe_energies == {}

    chem_utils.bst_for_each_contnr_no_opt([contnr], mols, 1, 1)

    # The three candidates have different formulas, so each is its own group,
    # and with thoroughness=1 each group has one candidate measured.
    assert len(contnr.probe_energies) == 3


class _ScoredConf:
    """A stand-in conformer carrying a fixed ranking energy."""

    def __init__(self, energy: float) -> None:
        self.energy = energy


def _scored_mol(smiles: str, energy: float) -> MyMol.MyMol:
    """Build a candidate whose ranking energy is fixed, so no embedding runs.

    first_conf_energy reads an existing conformer's energy directly, which
    lets a test choose the ranking without RDKit's embedder.

    Args:
        smiles: SMILES of the candidate.
        energy: The energy it should rank with.

    Returns:
        A MyMol carrying one stand-in conformer.
    """
    mol = MyMol.MyMol(smiles)
    mol.conformers = [_ScoredConf(energy)]
    return mol


def test_comparable_group_key_separates_protonation_states_only() -> None:
    """Only molecules with the same atoms may share a ranking group.

    Tautomers and stereoisomers keep their formula and charge, so they group
    with the protonation state they came from; different protonation states
    must not.
    """
    key = chem_utils._comparable_group_key

    assert key(MyMol.MyMol("CC(=O)O")) != key(MyMol.MyMol("CC(=O)[O-]"))
    assert key(MyMol.MyMol("CC(C)=O")) == key(MyMol.MyMol("C=C(C)O"))
    assert key(MyMol.MyMol("C[C@H](N)C(=O)O")) == key(MyMol.MyMol("C[C@@H](N)C(=O)O"))


def test_pick_lowest_enrgy_mols_keeps_every_protonation_state() -> None:
    """A higher-energy protonation state must still get a slot.

    Regression: candidates were ranked on one pooled list, so three neutral
    C2H4O2 isomers outranked the acetate anion, whose UFF energy is not on the
    same scale. Each group now gets a slot before any group gets a second.
    """
    neutral = [
        _scored_mol("CC(=O)O", 1.0),
        _scored_mol("OCC=O", 2.0),
        _scored_mol("COC=O", 3.0),
    ]
    anion = _scored_mol("CC(=O)[O-]", 100.0)

    for _ in range(20):
        kept = chem_utils.pick_lowest_enrgy_mols(neutral + [anion], 3, 3)
        assert {m.smiles() for m in kept} == {
            neutral[0].smiles(),
            neutral[1].smiles(),
            anion.smiles(),
        }


def test_pick_lowest_enrgy_mols_samples_every_group() -> None:
    """Sampling must not eliminate a protonation state before it is scored.

    Regression: the num * thoroughness sample was drawn from the pooled list,
    so a group with one member was often left out of the sample entirely. The
    sample is now drawn per group.
    """
    hexanes = [
        _scored_mol(smi, float(i))
        for i, smi in enumerate(
            ["CCCCCC", "CCCC(C)C", "CCC(C)CC", "CC(C)C(C)C", "CCC(C)(C)C"]
        )
    ]
    pentane = _scored_mol("CCCCC", 100.0)

    for _ in range(20):
        kept = chem_utils.pick_lowest_enrgy_mols(hexanes + [pentane], 2, 1)
        assert pentane.smiles() in {m.smiles() for m in kept}
