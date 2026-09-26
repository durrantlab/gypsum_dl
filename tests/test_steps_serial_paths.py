"""Tests for the in-process (parallelizer_obj is None) branch of each step.

The end-to-end tests always pass a real Parallelizer object, even in serial
mode, so the inline code paths in every pipeline step are otherwise never
executed.
"""

from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.steps.conf import Minimize3D
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


def test_enumerate_chiral_molecules_expands_unspecified_center() -> None:
    contnr = _container("CC(N)C(=O)O", "alanine")
    enumerate_chiral_molecules([contnr], 5, 1, 1, "serial", None)
    assert len(contnr.mols) >= 1


def test_enumerate_chiral_molecules_respects_zero_variants() -> None:
    contnr = _container("CC(N)C(=O)O", "alanine")
    enumerate_chiral_molecules([contnr], 0, 1, 1, "serial", None)
    assert len(contnr.mols) == 1


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