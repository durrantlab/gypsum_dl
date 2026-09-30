"""Tests for the in-process (parallelizer_obj is None) branch of each step.

The end-to-end tests always pass a real Parallelizer object, even in serial
mode, so the inline code paths in every pipeline step are otherwise never
executed.
"""

import numpy

from gypsum_dl import chem_utils
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.MyMol import ContnrProps, MyMol
from gypsum_dl.parallelizer import Parallelizer, seed_generators
from gypsum_dl.steps.conf import Minimize3D
from gypsum_dl.steps.smiles import (
    AddHydrogens,
    EnumerateChiralMols,
    EnumerateDoubleBonds,
    MakeTautomers,
)
from gypsum_dl.steps.conf.Convert2DTo3D import convert_2d_to_3d
from gypsum_dl.steps.conf.GenerateAlternate3DNonaromaticRingConfs import (
    generate_alternate_3d_nonaromatic_ring_confs,
    parallel_get_ring_confs,
)
from gypsum_dl.steps.conf.Minimize3D import minimize_3d
from gypsum_dl.steps.smiles.AddHydrogens import add_hydrogens
from gypsum_dl.steps.smiles.DeSaltOrigSmiles import desalt_orig_smi
from gypsum_dl.steps.smiles.DurrantLabFilter import (
    durrant_lab_contains_bad_substr,
    durrant_lab_filters,
)
from gypsum_dl.steps.smiles.EnumerateChiralMols import enumerate_chiral_molecules
from gypsum_dl.steps.smiles.EnumerateDoubleBonds import enumerate_double_bonds
from gypsum_dl.steps.smiles.MakeTautomers import make_tauts


def _container(smiles: str, name: str) -> MolContainer:
    """Build a single-molecule container ready for a pipeline step.

    Steps expect containers that already hold at least one variant, which is
    normally the job of the desalting step.

    Args:
        smiles: SMILES string for the input molecule.
        name: Ligand name.

    Returns:
        A populated MolContainer at container index zero.
    """
    contnr = MolContainer(smiles, name, 0, {})
    contnr.add_smiles(smiles)
    return contnr


def _container_at_idx(smiles: str, name: str, idx: int) -> MolContainer:
    """Build a populated container whose contnr_idx is not its list position.

    Every step regroups its work by contnr_idx, and several of them used that
    value to index the container list directly. Handing a step a single
    container that sits at position zero but carries a higher contnr_idx is
    what separates the two conventions.

    Args:
        smiles: SMILES string for the input molecule.
        name: Ligand name.
        idx: Container index to assign.

    Returns:
        A MolContainer holding a single variant.
    """
    contnr = MolContainer(smiles, name, idx, {})
    contnr.add_smiles(smiles)
    return contnr


def test_add_hydrogens_finds_failed_container_by_index(monkeypatch) -> None:
    # Regression: fnd_contnrs_not_represntd reports contnr_idx values, but the
    # carry-over loop used them to index the container list. A container whose
    # index is not its position was then read from the wrong slot, or raised
    # IndexError.
    contnr = _container_at_idx("CC(=O)O", "acetic_acid", 3)
    monkeypatch.setattr(AddHydrogens, "parallel_add_H", lambda *a: [])

    add_hydrogens([contnr], 6.4, 8.4, 1.0, 5, 1, 1, "serial", None)

    assert len(contnr.mols) == 1
    assert contnr.mols[0].genealogy[-1] == (
        "(WARNING: Gypsum-DL could not assign ionization states)"
    )


def test_add_hydrogens_fallback_does_not_alias_the_reference_mol(monkeypatch) -> None:
    # Regression: the carry-over added failed_contnr.mol_orig_frm_inp_smi
    # itself to the working set and overwrote its genealogy, so the container
    # lost its record of the input (including the entries the desalter had
    # stamped) and later steps rewrote the record in place.
    contnr = _container_at_idx("CC(=O)O", "acetic_acid", 3)
    reference = contnr.mol_orig_frm_inp_smi
    reference.genealogy = ["CC(=O)O (source)"]
    monkeypatch.setattr(AddHydrogens, "parallel_add_H", lambda *a: [])

    add_hydrogens([contnr], 6.4, 8.4, 1.0, 5, 1, 1, "serial", None)

    assert len(contnr.mols) == 1
    assert contnr.mols[0] is not reference
    assert reference.genealogy == ["CC(=O)O (source)"]


def test_parallel_add_h_does_not_mutate_shared_settings() -> None:
    # Regression: one protonation-settings dict is built per step and aliased
    # into every task tuple, and the worker wrote the current molecule's
    # SMILES into it. Correct only because the in-process path happens to be
    # sequential and the parallel paths pickle a copy per task.
    contnr = _container("CC(=O)O", "acetic_acid")
    settings = {
        "ph_min": 6.4,
        "ph_max": 8.4,
        "precision": 1.0,
        "max_variants": 5,
    }
    before = dict(settings)

    AddHydrogens.parallel_add_H(contnr, settings)

    assert settings == before


def test_add_hydrogens_fallback_does_not_invent_a_desalting_step(monkeypatch) -> None:
    # Regression (F3): the fallback replaced the genealogy with three lines it
    # synthesized, one of them a "(desalted)" entry, so a molecule that was
    # never desalted came out of the run asserting a step that never ran.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    desalt_orig_smi([contnr])
    monkeypatch.setattr(AddHydrogens, "parallel_add_H", lambda *a: [])

    add_hydrogens([contnr], 6.4, 8.4, 1.0, 5, 1, 1, "serial", None)

    assert len(contnr.mols) == 1
    assert contnr.mols[0].genealogy == [
        "CCO (source)",
        "(WARNING: Gypsum-DL could not assign ionization states)",
    ]


def test_add_hydrogens_fallback_keeps_the_input_smiles(monkeypatch) -> None:
    # Regression (F3): update_orig_smi has by then replaced orig_smi with the
    # desalted SMILES, so the synthesized "(source)" and "(desalted)" lines
    # were the same string and the salt form the user supplied was gone from
    # the record the README tells them to read.
    contnr = MolContainer("CCCCCCO.C", "salt", 0, {})
    desalt_orig_smi([contnr])
    monkeypatch.setattr(AddHydrogens, "parallel_add_H", lambda *a: [])

    add_hydrogens([contnr], 6.4, 8.4, 1.0, 5, 1, 1, "serial", None)

    assert len(contnr.mols) == 1
    genealogy = contnr.mols[0].genealogy
    assert len(genealogy) == 3
    assert genealogy[0] == "CCCCCCO.C (source)"
    assert genealogy[1].endswith(" (desalted)")
    assert genealogy[1] != genealogy[0]
    assert genealogy[2] == "(WARNING: Gypsum-DL could not assign ionization states)"


def test_parallel_add_h_accepts_a_single_smiles_string(monkeypatch) -> None:
    # Regression (F4): the return of protonate_smiles was fed straight to a
    # list comprehension, so a release that handed back one SMILES rather than
    # a one-element sequence would be walked character by character. The
    # single-atom fragments that happen to parse ("C", "O") survive
    # remove_bizarre_substruc and the Durrant filters, so the run would
    # quietly produce a library of methane.
    contnr = _container("CC(=O)O", "acetic_acid")
    monkeypatch.setattr(AddHydrogens, "protonate_smiles", lambda **kwargs: "CC(=O)[O-]")
    settings = {
        "ph_min": 6.4,
        "ph_max": 8.4,
        "precision": 1.0,
        "max_variants": 5,
    }

    results = AddHydrogens.parallel_add_H(contnr, settings)

    assert len(results) == 1
    assert "[O-]" in results[0].smiles()


def test_parallel_add_h_actually_ionizes() -> None:
    # The existing coverage asserts only that a container ends up with at
    # least one variant, which the F3 fallback satisfies too, so it cannot
    # tell working ionization apart from none at all. A carboxylic acid well
    # above its pKa has to come back deprotonated.
    contnr = _container("CC(=O)O", "acetic_acid")
    settings = {
        "ph_min": 12.0,
        "ph_max": 12.0,
        "precision": 1.0,
        "max_variants": 5,
    }

    results = AddHydrogens.parallel_add_H(contnr, settings)

    assert any("[O-]" in mol.smiles() for mol in results)


def test_enumerate_chiral_finds_failed_container_by_index(monkeypatch) -> None:
    # Same index-versus-position confusion in the enantiomer carry-over.
    contnr = _container_at_idx("CC(N)C(=O)O", "alanine", 3)
    original_mol = contnr.mols[0]
    monkeypatch.setattr(EnumerateChiralMols, "parallel_get_chiral", lambda *a: None)

    enumerate_chiral_molecules([contnr], 5, 1, 1, "serial", None)

    assert original_mol.genealogy[-1] == "(WARNING: Unable to generate enantiomers)"
    assert contnr.mols == [original_mol]


def test_enumerate_double_bonds_finds_failed_container_by_index(monkeypatch) -> None:
    # Same index-versus-position confusion in the double-bond carry-over.
    contnr = _container_at_idx("CC=CCC", "pentene", 3)
    original_mol = contnr.mols[0]
    monkeypatch.setattr(
        EnumerateDoubleBonds, "parallel_get_double_bonded", lambda *a: None
    )

    enumerate_double_bonds([contnr], 5, 1, 1, "serial", None)

    assert (
        original_mol.genealogy[-1]
        == "(WARNING: Unable to generate double-bond variant)"
    )
    assert contnr.mols == [original_mol]


def test_minimize_3d_populates_container_by_index() -> None:
    # Regression: minimize_3d emptied and repopulated containers with
    # contnrs[mol.contnr_idx], which is only the right container while every
    # index matches its list position.
    contnr = _container_at_idx("CCO", "ethanol", 3)
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)

    minimize_3d([contnr], 1, 1, 1, False, "serial", None)

    assert len(contnr.mols) == 1
    assert "Energy" in contnr.mols[0].mol_props


def test_generate_alternate_ring_confs_by_index() -> None:
    # Regression: this step tracked ring-bearing containers by list position
    # but grouped its results by contnr_idx, then indexed the container list
    # with the grouped keys. All three had to agree.
    contnr = _container_at_idx("C1CCCCC1", "cyclohexane", 3)
    convert_2d_to_3d([contnr], 3, 2, 1, "serial", None)

    generate_alternate_3d_nonaromatic_ring_confs(
        [contnr], 3, 2, 1, False, "serial", None
    )

    assert len(contnr.mols) >= 1
    assert len(contnr.mols[0].conformers) == 1


def test_generate_alternate_ring_confs_flags_failure_by_index(monkeypatch) -> None:
    # The no-results branch subtracted the grouped keys (contnr_idx values)
    # from a set of list positions, so with the two conventions disagreeing it
    # flagged the wrong container or raised IndexError.
    contnr = _container_at_idx("C1CCCCC1", "cyclohexane", 3)
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    original_mol = contnr.mols[0]
    monkeypatch.setattr(
        "gypsum_dl.steps.conf.GenerateAlternate3DNonaromaticRingConfs."
        "parallel_get_ring_confs",
        lambda *a: None,
    )

    generate_alternate_3d_nonaromatic_ring_confs(
        [contnr], 1, 1, 1, False, "serial", None
    )

    assert original_mol.genealogy[-1] == (
        "(WARNING: Could not generate alternate conformations of nonaromatic ring)"
    )


def test_desalt_orig_smi_keeps_largest_fragment() -> None:
    contnr = MolContainer("CCCCCCO.C", "salt", 0, {})
    desalt_orig_smi([contnr])
    assert "." not in contnr.orig_smi
    assert len(contnr.mols) == 1


def test_desalt_orig_smi_leaves_single_fragment_alone() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    desalt_orig_smi([contnr])
    assert contnr.orig_smi == "CCO"
    assert len(contnr.mols) == 1


def test_add_hydrogens_generates_ionization_states() -> None:
    contnr = _container("CC(=O)O", "acetic_acid")
    add_hydrogens([contnr], 6.4, 8.4, 1.0, 5, 1, 1, "serial", None)
    assert len(contnr.mols) >= 1


def test_make_tauts_keeps_at_least_one_variant() -> None:
    contnr = _container("CC(=O)CC", "butanone")
    make_tauts([contnr], 5, 1, 1, "serial", False, None)
    assert len(contnr.mols) >= 1


def test_make_tauts_allows_chirality_changes_when_requested() -> None:
    contnr = _container("CC(=O)CC", "butanone")
    make_tauts([contnr], 5, 1, 1, "serial", True, None)
    assert len(contnr.mols) >= 1


def test_make_tauts_respects_zero_variants() -> None:
    contnr = _container("CC(=O)CC", "butanone")
    make_tauts([contnr], 0, 1, 1, "serial", False, None)
    assert len(contnr.mols) == 1


def test_make_tauts_scales_enumerator_budget_by_thoroughness(monkeypatch) -> None:
    # Regression: make_tauts handed MolVS the bare variant cap, so with the
    # defaults the enumerator stopped expanding at 5 while every other step
    # generated thoroughness * max_variants_per_compound candidates. MolVS uses
    # max_tautomers as a stopping size for its breadth-first expansion, so the
    # bare cap also meant max_variants_per_compound=1 disabled tautomerization
    # outright (its loop body never runs at 1).
    budgets: list[int] = []

    class _RecordingEnumerator:
        """Capture the cap MolVS is constructed with, then stand in for it.

        The real enumerator's output is irrelevant here; what matters is the
        number make_tauts chose, which is otherwise invisible from outside.

        Args:
            max_tautomers: The stopping size under test.
        """

        def __init__(self, max_tautomers: int) -> None:
            budgets.append(max_tautomers)

        def enumerate(self, mol: object) -> list[object]:
            """Return the input unchanged, standing in for a real enumeration.

            Args:
                mol: The kekulized RDKit molecule handed to MolVS.

            Returns:
                A single-member list holding that same molecule.
            """
            return [mol]

    monkeypatch.setattr(
        MakeTautomers.tautomer, "TautomerEnumerator", _RecordingEnumerator
    )

    contnr = _container("CC(=O)CC", "butanone")
    make_tauts([contnr], 5, 3, 1, "serial", False, None)
    assert budgets == [15]

    budgets.clear()
    contnr = _container("CC(=O)CC", "butanone")
    make_tauts([contnr], 1, 3, 1, "serial", False, None)
    assert budgets == [3]


_REAL_MAKE_TAUT = MakeTautomers.parallel_make_taut


def _make_taut_failing_on_ethanol(
    mol: MyMol, props: ContnrProps, max_tauts: int
) -> list[MyMol] | None:
    """Stand in for parallel_make_taut, raising for one chosen compound.

    The failure has to originate inside the worker function itself so that it
    travels whichever dispatch path job_manager selects. The real function is
    bound at import time so the monkeypatched name cannot recurse into itself,
    and this lives at module scope so worker processes can unpickle it.

    Args:
        mol: The variant being tautomerized.
        props: The container-level fields describing the input compound.
        max_tauts: Size at which MolVS stops expanding the tautomer set.

    Returns:
        Whatever the real function returns, for every compound but ethanol.

    Raises:
        RuntimeError: For the compound named "ethanol".
    """
    if props["name"] == "ethanol":
        raise RuntimeError("simulated tautomerization failure")
    return _REAL_MAKE_TAUT(mol, props, max_tauts)


def _taut_test_contnrs() -> list[MolContainer]:
    """Build two distinctly indexed containers for the taut-failure tests.

    Every step regroups its results by contnr_idx, so a dropped molecule is
    only attributable to a container when the indices differ.

    Returns:
        Containers for ethanol (index 0) and butanone (index 1), each holding
        a single variant.
    """
    specs = [("CCO", "ethanol"), ("CC(=O)CC", "butanone")]
    contnrs = []
    for idx, (smiles, name) in enumerate(specs):
        contnr = MolContainer(smiles, name, idx, {})
        contnr.add_smiles(smiles)
        contnrs.append(contnr)
    return contnrs


def _contnr_smiles(contnrs: list[MolContainer]) -> list[list[str]]:
    """Summarize container membership so two runs can be compared directly.

    Args:
        contnrs: The containers to summarize.

    Returns:
        One sorted list of canonical SMILES per container.
    """
    return [sorted(mol.smiles() for mol in contnr.mols) for contnr in contnrs]


def test_make_tauts_drops_raising_molecule_in_process(monkeypatch) -> None:
    # Regression: the in-process branch called the worker function bare, so a
    # single molecule that raised aborted the whole run instead of being
    # dropped and carried over from the previous step.
    monkeypatch.setattr(
        MakeTautomers, "parallel_make_taut", _make_taut_failing_on_ethanol
    )
    contnrs = _taut_test_contnrs()
    before = _contnr_smiles(contnrs)

    make_tauts(contnrs, 50, 1, 1, "serial", False, None)

    assert _contnr_smiles(contnrs)[0] == before[0]
    assert len(contnrs[1].mols) >= 1


def test_make_tauts_failure_handling_matches_across_job_managers(monkeypatch) -> None:
    # Regression: the multiprocessing worker reported a raised exception as a
    # None result while both serial paths let it propagate, so the output of a
    # run depended on the job_manager setting.
    monkeypatch.setattr(
        MakeTautomers, "parallel_make_taut", _make_taut_failing_on_ethanol
    )

    in_process = _taut_test_contnrs()
    make_tauts(in_process, 50, 1, 1, "serial", False, None)

    serial_par = Parallelizer("serial", 1)
    serial = _taut_test_contnrs()
    make_tauts(serial, 50, 1, 1, "serial", False, serial_par)
    serial_par.end()

    mp_par = Parallelizer("multiprocessing", 2, True)
    multiproc = _taut_test_contnrs()
    make_tauts(multiproc, 50, 1, 2, "multiprocessing", False, mp_par)
    mp_par.end()

    assert _contnr_smiles(in_process) == _contnr_smiles(serial)
    assert _contnr_smiles(serial) == _contnr_smiles(multiproc)


def test_enumerate_chiral_molecules_expands_unspecified_center() -> None:
    contnr = _container("CC(N)C(=O)O", "alanine")
    enumerate_chiral_molecules([contnr], 5, 1, 1, "serial", None)
    assert len(contnr.mols) >= 1


def test_enumerate_chiral_molecules_respects_zero_variants() -> None:
    contnr = _container("CC(N)C(=O)O", "alanine")
    enumerate_chiral_molecules([contnr], 0, 1, 1, "serial", None)
    assert len(contnr.mols) == 1


def _capture_carried_over(monkeypatch, module):
    """Capture the flat list handed to bst_for_each_contnr_no_opt.

    Also stops the step from touching container membership so the assertion
    isn't masked by the crry_ovr_frm_lst_step_if_no_fnd default, which would
    otherwise repopulate the container regardless of which list got the mol.
    """
    captured = {}

    def fake_bst(contnrs, mol_lst, *args, **kwargs):
        captured["flat"] = mol_lst

    monkeypatch.setattr(chem_utils, "bst_for_each_contnr_no_opt", fake_bst)
    return captured


def test_enumerate_chiral_carries_over_failed_container(monkeypatch) -> None:
    # Regression: when a container yields no enantiomers, its existing mol must
    # be put back into the list passed downstream (flat), not the discarded one.
    contnr = _container("CC(N)C(=O)O", "alanine")
    original_mol = contnr.mols[0]
    monkeypatch.setattr(EnumerateChiralMols, "parallel_get_chiral", lambda *a: None)
    captured = _capture_carried_over(monkeypatch, EnumerateChiralMols)

    enumerate_chiral_molecules([contnr], 5, 1, 1, "serial", None)

    assert original_mol in captured["flat"]
    assert original_mol.genealogy[-1] == "(WARNING: Unable to generate enantiomers)"


def test_enumerate_double_bonds_carries_over_failed_container(monkeypatch) -> None:
    # Regression: same carry-over path for double-bond enumeration.
    contnr = _container("CC=CCC", "pentene")
    original_mol = contnr.mols[0]
    monkeypatch.setattr(
        EnumerateDoubleBonds, "parallel_get_double_bonded", lambda *a: None
    )
    captured = _capture_carried_over(monkeypatch, EnumerateDoubleBonds)

    enumerate_double_bonds([contnr], 5, 1, 1, "serial", None)

    assert original_mol in captured["flat"]
    assert (
        original_mol.genealogy[-1]
        == "(WARNING: Unable to generate double-bond variant)"
    )


def test_enumerate_double_bonds_expands_unspecified_bond() -> None:
    contnr = _container("CC=CCC", "pentene")
    enumerate_double_bonds([contnr], 5, 1, 1, "serial", None)
    assert len(contnr.mols) >= 1


def test_enumerate_double_bonds_enumerates_each_repeated_input() -> None:
    # Regression: the step deduplicated its variants over every container at
    # once, so when two inputs shared a SMILES (a repeated entry, or a free
    # acid and its salt after desalting) the first container kept the variants
    # and the second got none. It then fell back to its unenumerated mol with
    # a warning about conformers instead.
    first = _container_at_idx("CC=CC", "butene", 0)
    second = _container_at_idx("CC=CC", "butene_again", 1)

    enumerate_double_bonds([first, second], 2, 1, 1, "serial", None)

    for contnr in (first, second):
        smis = [m.smiles(True) for m in contnr.mols]
        assert len(smis) == 2, f"{contnr.name} ended with {smis}"
        assert all("/" in s or "\\" in s for s in smis)


def test_enumerate_double_bonds_respects_zero_variants() -> None:
    contnr = _container("CC=CCC", "pentene")
    enumerate_double_bonds([contnr], 0, 1, 1, "serial", None)
    assert len(contnr.mols) == 1


def test_durrant_lab_contains_bad_substr_detects_metals() -> None:
    assert durrant_lab_contains_bad_substr("CC(=O)[O-][Zn+2]") is True
    assert durrant_lab_contains_bad_substr("CCO") is False


def test_durrant_lab_filters_discards_boron() -> None:
    contnr = _container("B(O)(O)O", "boric_acid")
    durrant_lab_filters([contnr], 1, "serial", None)
    assert contnr.mols == []


def test_durrant_lab_filters_discards_metal_by_substring() -> None:
    # The metal check is a SMILES substring test rather than a substructure
    # match, so it is the one branch of the filter that no pattern covers.
    contnr = _container("CCO.[Zn+2]", "zinc_salt")
    durrant_lab_filters([contnr], 1, "serial", None)
    assert contnr.mols == []


def test_durrant_lab_filters_keeps_clean_molecule() -> None:
    contnr = _container("CCO", "ethanol")
    durrant_lab_filters([contnr], 1, "serial", None)
    assert len(contnr.mols) == 1


def test_durrant_lab_filters_sizes_its_variant_cap_to_the_candidates(
    monkeypatch,
) -> None:
    # Regression: the step passed a fixed cap of 1000 to the energy-based
    # selector, which is a no-op only while every container holds fewer
    # variants than that. max_variants_per_compound has no upper bound, so a
    # container allowed past the constant would have been silently pruned by a
    # step that is only supposed to rebuild container membership, with
    # throwaway conformers embedded to rank the variants it dropped. The cap is
    # now sized to the candidate list, so no container can exceed it.
    first = MolContainer("CCO", "ethanol", 0, {})
    first.add_smiles(["CCO", "CCCO", "CCCCO"])
    second = MolContainer("CCC", "propane", 1, {})
    second.add_smiles(["CCC", "CCCC"])
    candidates = len(first.mols) + len(second.mols)

    caps: list[int] = []
    real_pick = chem_utils.pick_lowest_enrgy_mols

    def recording_pick(
        mol_lst: list[MyMol],
        num: int,
        thoroughness: int,
        energy_cache: dict[str, float | None] | None = None,
    ) -> list[MyMol]:
        """Record the cap the filter asks for, then select as usual.

        Args:
            mol_lst: The candidate variants for one container.
            num: The number of variants the caller is willing to keep.
            thoroughness: How many candidates to evaluate per kept variant.
            energy_cache: The container's ranking-energy cache, passed through
                untouched.

        Returns:
            The variants the real selector would keep.
        """
        caps.append(num)
        return real_pick(mol_lst, num, thoroughness, energy_cache)

    monkeypatch.setattr(chem_utils, "pick_lowest_enrgy_mols", recording_pick)

    durrant_lab_filters([first, second], 1, "serial", None)

    assert caps == [candidates, candidates]
    assert len(first.mols) == 3
    assert len(second.mols) == 2
    # Nothing was ranked, so nothing was embedded.
    assert all(not mol.conformers for contnr in (first, second) for mol in contnr.mols)


def test_durrant_lab_filters_does_not_claim_originals_were_kept(capsys) -> None:
    # Regression (B14): the step empties every container before calling
    # bst_for_each_contnr_no_opt, so the carry-over default logged "Keeping
    # original conformers" for a container it had just emptied, contradicting
    # the per-variant discard message and gypsum_dl_failed.smi.
    contnr = _container("B(O)(O)O", "boric_acid")
    durrant_lab_filters([contnr], 1, "serial", None)
    log = " ".join(capsys.readouterr().out.split())
    assert contnr.mols == []
    assert "discarding it" in log
    assert "Keeping original" not in log


def test_durrant_lab_filters_discards_variant_with_no_rdkit_mol() -> None:
    # Regression (B14): HasSubstructMatch was called without a None guard, so a
    # variant whose RDKit mol failed to build raised AttributeError (a silent
    # drop under multiprocessing).
    contnr = _container("CCO", "ethanol")
    contnr.mols[0].rdkit_mol = None
    durrant_lab_filters([contnr], 1, "serial", None)
    assert contnr.mols == []


def test_durrant_lab_filters_keeps_siblings_of_a_none_rdkit_mol() -> None:
    contnr = _container("CCO", "ethanol")
    contnr.add_smiles("CCCO")
    assert len(contnr.mols) == 2
    contnr.mols[1].rdkit_mol = None
    durrant_lab_filters([contnr], 1, "serial", None)
    assert [m.smiles() for m in contnr.mols] == ["CCO"]


def test_convert_2d_to_3d_assigns_a_conformer() -> None:
    contnr = _container("CCO", "ethanol")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    assert len(contnr.mols) == 1
    assert len(contnr.mols[0].conformers) == 1


def test_convert_2d_to_3d_discards_bizarre_substructure() -> None:
    contnr = _container("CC[CH2-]", "carbanion")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    assert contnr.mols == []


def _seeded_conversion_coords(parallelizer_obj: Parallelizer | None) -> list:
    """Convert one molecule to 3D under a fixed seed and report the geometry.

    The seed is installed immediately before the step so that both dispatch
    paths start from the same generator state; anything that reached a
    different state by the time the embedding drew its RDKit seed would show
    up as different coordinates.

    Args:
        parallelizer_obj: The Parallelizer to dispatch through, or None to
            exercise the in-process branch.

    Returns:
        The first conformer's coordinates, as nested lists so two runs can be
            compared directly.
    """
    seed_generators(1234)
    contnr = _container("CCCCO", "butanol")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", parallelizer_obj)
    return contnr.mols[0].conformers[0].coords().tolist()


def test_in_process_dispatch_is_seeded_like_serial_mode() -> None:
    # Regression: the in-process branch of every step called run_one with no
    # seed, so --random_seed did not reach it. That branch is what an mpi
    # worker and any library caller passing parallelizer_obj=None run, which
    # left those runs unseeded while serial and multiprocessing runs were
    # reproducible.
    in_process = _seeded_conversion_coords(None)

    serial_par = Parallelizer("serial", 1)
    serial = _seeded_conversion_coords(serial_par)
    serial_par.end()

    assert in_process == serial


def test_in_process_dispatch_repeats_itself_under_one_seed() -> None:
    # Control for the test above: the agreement is between two reproducible
    # runs, not between two runs that both happen to be unseeded.
    assert _seeded_conversion_coords(None) == _seeded_conversion_coords(None)


def test_generate_alternate_ring_confs_keeps_variants() -> None:
    contnr = _container("C1CCCCC1", "cyclohexane")
    convert_2d_to_3d([contnr], 3, 2, 1, "serial", None)
    generate_alternate_3d_nonaromatic_ring_confs(
        [contnr], 3, 2, 1, False, "serial", None
    )
    assert len(contnr.mols) >= 1
    assert len(contnr.mols[0].conformers) == 1


class _TiedConformer:
    """A stand-in conformer that only has to report an energy.

    The tie-breaking code reads `conformers[0].energy` and nothing else, and
    the tie that matters in practice is the sentinel energy MyConformer
    assigns when the UFF setup fails, which no real conformer can be made to
    produce on demand.
    """

    def __init__(self, energy: float) -> None:
        self.energy = energy


def _tied_ring_conf_variant(smiles: str) -> MyMol:
    """Build a ring-conformer result whose energy ties with its siblings.

    Args:
        smiles: SMILES string for the variant.

    Returns:
        A MyMol at container index zero carrying the sentinel energy.
    """
    mol = MyMol(smiles)
    mol.contnr_idx = 0
    mol.conformers = [_TiedConformer(9999)]
    return mol


def test_generate_alternate_ring_confs_breaks_energy_ties_by_smiles(
    monkeypatch,
) -> None:
    # Regression: the (energy, mol) pairs were sorted with a bare sort(), so
    # tied energies fell through to comparing the MyMol objects themselves,
    # which compare by hash(canonical_smiles). CPython salts string hashing
    # per invocation, and ties are routine (a failed UFF setup gives every
    # conformer the same sentinel energy), so which variants survived the
    # max_variants_per_compound trim changed between otherwise identical runs,
    # including runs with random_seed set. The expected order below is the
    # SMILES order, which does not move with the hash seed.
    contnr = _container("C1CCCCC1", "cyclohexane")
    variants = [_tied_ring_conf_variant(smi) for smi in ("CCCO", "CCO", "CCCCO")]
    monkeypatch.setattr(
        "gypsum_dl.steps.conf.GenerateAlternate3DNonaromaticRingConfs."
        "parallel_get_ring_confs",
        lambda *a: variants,
    )

    generate_alternate_3d_nonaromatic_ring_confs(
        [contnr], 2, 1, 1, False, "serial", None
    )

    assert [mol.smiles() for mol in contnr.mols] == ["CCCCO", "CCCO"]


def test_generate_alternate_ring_confs_tolerates_a_tie_with_no_smiles(
    monkeypatch,
) -> None:
    # smiles() reports failure as None, which the tie-break key has to absorb:
    # comparing None against a string raises TypeError, and this sort runs in
    # the main process after the step has fanned out.
    contnr = _container("C1CCCCC1", "cyclohexane")
    variants = [_tied_ring_conf_variant(smi) for smi in ("CCO", "CCCO")]
    variants[0].can_smi = None
    monkeypatch.setattr(
        "gypsum_dl.steps.conf.GenerateAlternate3DNonaromaticRingConfs."
        "parallel_get_ring_confs",
        lambda *a: variants,
    )

    generate_alternate_3d_nonaromatic_ring_confs(
        [contnr], 2, 1, 1, False, "serial", None
    )

    assert len(contnr.mols) == 2


def test_parallel_get_ring_confs_clusters_several_ring_conformers() -> None:
    # Regression: the degenerate-shape guard tested pts.shape == (0,), a shape
    # the earlier "no rings" return already makes unreachable, so it was
    # protecting nothing. Written as a check on pts.shape[0] instead it would
    # have swallowed every molecule, since the rmsd points are measured
    # relative to the first conformer and so are one row short until the
    # vstack adds it. Cyclohexane keeps several distinct ring geometries, so
    # collapsing to a single conformer here is a real loss of output.
    contnr = _container("C1CCCCC1", "cyclohexane")
    convert_2d_to_3d([contnr], 3, 2, 1, "serial", None)
    mol = contnr.mols[0]

    results = parallel_get_ring_confs(mol, 3, 2, False)

    # Precondition: with only one surviving conformer a single variant would
    # be the correct answer, and the assertion below would say nothing.
    assert len(mol.conformers) > 1
    assert len(results) > 1
    coords = {
        tuple(numpy.round(variant.conformers[0].coords().flatten(), 3))
        for variant in results
    }
    assert len(coords) == len(results)


def test_parallel_get_ring_confs_keeps_a_lone_conformer() -> None:
    # A single conformer leaves every ring rmsd list empty, so pts arrives
    # with shape (0, num_rings). That is the reference conformer rather than a
    # molecule with nothing to cluster, and the guard has to let it reach the
    # vstack that supplies its row of zeros.
    contnr = _container("C1CCCCC1", "cyclohexane")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    mol = contnr.mols[0]

    results = parallel_get_ring_confs(mol, 1, 1, False)

    assert len(mol.conformers) == 1
    assert len(results) == 1
    assert len(results[0].conformers) == 1


def test_generate_alternate_ring_confs_honors_the_minimize_flag(monkeypatch) -> None:
    # Regression: this step minimizes every conformer it generates, with a
    # hardcoded True, so --skip_optimize_geometry silently did not apply to any
    # molecule with a non-aromatic ring. The flag has to reach add_conformers,
    # which means travelling through the parallelizer params tuple as well.
    contnr = _container("C1CCCCC1", "cyclohexane")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)

    captured: list[bool] = []
    original_add_conformers = MyMol.add_conformers

    def spy(
        self: MyMol, num: int, rmsd_cutoff: float = 0.1, minimize: bool = True
    ) -> None:
        """Record the minimize argument, then do the real work.

        The generated conformers are rebuilt into fresh MyConformer objects
        before the step returns, and those always report minimized == False, so
        the output cannot show whether minimization happened.

        Args:
            self: The molecule gaining conformers.
            num: Number of conformers requested.
            rmsd_cutoff: Redundancy cutoff, as in MyMol.add_conformers.
            minimize: Whether to minimize the conformers.
        """
        captured.append(minimize)
        original_add_conformers(self, num, rmsd_cutoff, minimize)

    monkeypatch.setattr(MyMol, "add_conformers", spy)

    generate_alternate_3d_nonaromatic_ring_confs(
        [contnr], 1, 1, 1, False, "serial", None, minimize=False
    )

    assert captured == [False]


def test_generate_alternate_ring_confs_skips_molecules_without_rings() -> None:
    contnr = _container("CCO", "ethanol")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    generate_alternate_3d_nonaromatic_ring_confs(
        [contnr], 1, 1, 1, False, "serial", None
    )
    assert len(contnr.mols) == 1


def test_minimize_3d_records_energy() -> None:
    contnr = _container("CCO", "ethanol")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    minimize_3d([contnr], 1, 1, 1, False, "serial", None)
    assert "Energy" in contnr.mols[0].mol_props


def test_minimize_3d_survives_none_worker_result(monkeypatch) -> None:
    # A worker returns None for a molecule with no acceptable conformers.
    # minimize_3d must skip it rather than dereference None, and still place
    # the surviving molecule in its container. Distinct indices because the
    # stub below keys on contnr_idx: built with two containers at index zero,
    # nothing was ever skipped and both minimized mols landed in the first
    # container.
    contnrs = [
        _container("CCO", "ethanol"),
        _container_at_idx("CCC", "propane", 1),
    ]
    for contnr in contnrs:
        convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)

    real_parallel_minit = Minimize3D.parallel_minit

    def stub(mol, *args):
        if mol.contnr_idx == 1:
            return None
        return real_parallel_minit(mol, *args)

    monkeypatch.setattr(Minimize3D, "parallel_minit", stub)

    minimize_3d(contnrs, 1, 1, 1, False, "serial", None)

    # The surviving molecule is minimized and recorded; the None result is
    # skipped without raising, leaving its container's pre-min mol untouched.
    assert "Energy" in contnrs[0].mols[0].mol_props
    assert len(contnrs[1].mols) == 1


def test_minimize_3d_skips_ring_mols_by_default() -> None:
    # Containers with non-aromatic rings are minimized by the ring-conformer
    # step, so minimize_3d leaves them alone unless asked.
    contnr = _container("C1CCCCC1", "cyclohexane")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    minimize_3d([contnr], 1, 1, 1, False, "serial", None)
    assert "Energy" not in contnr.mols[0].mol_props


def test_minimize_3d_records_energy_for_ring_mols_when_requested() -> None:
    # Regression: with the ring-conformer step skipped, nothing else minimizes
    # molecules with non-aromatic rings, so minimize_3d has to.
    contnr = _container("C1CCCCC1", "cyclohexane")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    minimize_3d([contnr], 1, 1, 1, False, "serial", None, include_nonaro_rings=True)
    assert len(contnr.mols) == 1
    assert "Energy" in contnr.mols[0].mol_props
    assert len(contnr.mols[0].conformers) == 1
