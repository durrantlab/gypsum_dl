"""Unit tests for MolContainer bookkeeping."""

import pytest

from gypsum_dl import MyMol
from gypsum_dl.MolContainer import MolContainer


def test_container_records_structural_counts() -> None:
    contnr = MolContainer("C[C@H](N)C(=O)O", "alanine", 0, {"src": "test"})
    assert contnr.orig_smi_canonical
    assert contnr.num_specif_chiral_cntrs == 1
    assert contnr.num_unspecif_chiral_cntrs == 1
    assert contnr.num_nonaro_rngs == 0
    assert contnr.carbon_hydrogen_count == 4


def test_container_counts_nonaromatic_rings() -> None:
    assert MolContainer("C1CCCCC1", "cyclohexane", 0, {}).num_nonaro_rngs == 1


def test_add_smiles_skips_duplicates() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles(["CCO", "OCC", "CCC"])
    assert len(contnr.mols) == 2


def test_add_smiles_preserves_name_and_orig_smi() -> None:
    # Regression: add_smiles must copy the container name and orig_smi onto the
    # variant, not overwrite the name with orig_smi.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("OCC")
    mol = contnr.mols[0]
    assert mol.name == "ethanol"
    assert mol.orig_smi == "CCO"


def test_mol_with_smiles_is_in_contnr_detects_existing() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    assert contnr.mol_with_smiles_is_in_contnr("OCC") is True


def test_add_container_properties_copies_properties() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {"activity": "1.0"})
    contnr.add_smiles("CCO")
    contnr.add_container_properties()
    assert contnr.mols[0].mol_props["activity"] == "1.0"
    assert contnr.mols[0].rdkit_mol.GetProp("activity") == "1.0"


def test_remove_identical_mols_from_contnr() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_mol(MyMol.MyMol("CCO"))
    contnr.add_mol(MyMol.MyMol("OCC"))
    contnr.remove_identical_mols_from_contnr()
    assert len(contnr.mols) == 1


def test_all_can_noh_smiles() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    assert contnr.all_can_noh_smiles() == ["CCO"]


def test_get_frags_of_orig_smi_caches_result() -> None:
    contnr = MolContainer("CCO.CC", "salt", 0, {})
    frags = contnr.get_frags_of_orig_smi()
    assert len(frags) == 2
    assert contnr.get_frags_of_orig_smi() is frags


def test_update_orig_smi_resets_state() -> None:
    contnr = MolContainer("CCO.CC", "salt", 0, {})
    contnr.add_smiles("CCO")
    contnr.update_orig_smi("CCO")
    assert contnr.mols == []
    assert contnr.orig_smi == "CCO"
    assert contnr.orig_smi_deslt == "CCO"
    assert len(contnr.get_frags_of_orig_smi()) == 1


def test_update_idx_propagates_to_original_mol() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.update_idx(5)
    assert contnr.contnr_idx == 5
    assert contnr.mol_orig_frm_inp_smi.contnr_idx == 5


def test_update_idx_requires_int() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    with pytest.raises(Exception):
        contnr.update_idx("5")