"""Tests for the in-process (parallelizer_obj is None) branch of each step.

The end-to-end tests always pass a real Parallelizer object, even in serial
mode, so the inline code paths in every pipeline step are otherwise never
executed.
"""

from gypsum_dl import chem_utils
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.MyMol import MyMol
from gypsum_dl.parallelizer import Parallelizer
from gypsum_dl.steps.conf import Minimize3D
from gypsum_dl.steps.smiles import (
    EnumerateChiralMols,
    EnumerateDoubleBonds,
    MakeTautomers,
)
from gypsum_dl.steps.conf.Convert2DTo3D import convert_2d_to_3d
from gypsum_dl.steps.conf.GenerateAlternate3DNonaromaticRingConfs import (
    generate_alternate_3d_nonaromatic_ring_confs,
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


def test_desalt_orig_smi_keeps_largest_fragment() -> None:
    contnr = MolContainer("CCCCCCO.C", "salt", 0, {})
    desalt_orig_smi([contnr], 1, "serial", None)
    assert "." not in contnr.orig_smi
    assert len(contnr.mols) == 1


def test_desalt_orig_smi_leaves_single_fragment_alone() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    desalt_orig_smi([contnr], 1, "serial", None)
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


_REAL_MAKE_TAUT = MakeTautomers.parallel_make_taut


def _make_taut_failing_on_ethanol(
    contnr: MolContainer, mol_index: int, max_variants_per_compound: int
) -> list[MyMol] | None:
    """Stand in for parallel_make_taut, raising for one chosen container.

    The failure has to originate inside the worker function itself so that it
    travels whichever dispatch path job_manager selects. The real function is
    bound at import time so the monkeypatched name cannot recurse into itself,
    and this lives at module scope so worker processes can unpickle it.

    Args:
        contnr: The molecule container being tautomerized.
        mol_index: Index of the molecule within the container.
        max_variants_per_compound: Cap on the number of tautomers enumerated.

    Returns:
        Whatever the real function returns, for every container but ethanol.

    Raises:
        RuntimeError: For the container named "ethanol".
    """
    if contnr.name == "ethanol":
        raise RuntimeError("simulated tautomerization failure")
    return _REAL_MAKE_TAUT(contnr, mol_index, max_variants_per_compound)


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
    monkeypatch.setattr(EnumerateDoubleBonds, "parallel_get_double_bonded", lambda *a: None)
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


def test_durrant_lab_filters_keeps_clean_molecule() -> None:
    contnr = _container("CCO", "ethanol")
    durrant_lab_filters([contnr], 1, "serial", None)
    assert len(contnr.mols) == 1


def test_convert_2d_to_3d_assigns_a_conformer() -> None:
    contnr = _container("CCO", "ethanol")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    assert len(contnr.mols) == 1
    assert len(contnr.mols[0].conformers) == 1


def test_convert_2d_to_3d_discards_bizarre_substructure() -> None:
    contnr = _container("CC[CH2-]", "carbanion")
    convert_2d_to_3d([contnr], 1, 1, 1, "serial", None)
    assert contnr.mols == []


def test_generate_alternate_ring_confs_keeps_variants() -> None:
    contnr = _container("C1CCCCC1", "cyclohexane")
    convert_2d_to_3d([contnr], 3, 2, 1, "serial", None)
    generate_alternate_3d_nonaromatic_ring_confs(
        [contnr], 3, 2, 1, False, "serial", None
    )
    assert len(contnr.mols) >= 1
    assert len(contnr.mols[0].conformers) == 1


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
    # the surviving molecule in its container.
    contnrs = [_container("CCO", "ethanol"), _container("CCC", "propane")]
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
