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


def test_mol_with_smiles_is_in_contnr_is_a_predicate() -> None:
    # Regression: the method returned either True or a freshly built MyMol, so
    # add_smiles decided membership with `result != True`, which ran
    # MyMol.__eq__ against a bool and compared canonical-SMILES hashes with
    # hash(True) == 1. Both answers must now be plain booleans.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")

    assert contnr.mol_with_smiles_is_in_contnr("CCC") is False
    assert contnr.mol_with_smiles_is_in_contnr("OCC") is True


def test_contains_canonical_smiles_treats_unknown_as_absent() -> None:
    # Two molecules that both failed to canonicalize have not been shown to be
    # the same molecule, so a non-string canonical SMILES must never match an
    # existing entry (which would silently drop the new variant).
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    contnr.mols[0].can_smi = None

    assert contnr.contains_canonical_smiles(None) is False
    assert contnr.contains_canonical_smiles("CCO") is False


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


def test_update_orig_smi_derives_the_same_fields_as_construction() -> None:
    # Regression: the constructor and update_orig_smi each carried their own
    # copy of the same derivation, and they drifted (one refreshed a derived
    # field that the other left describing the pre-desalt molecule). Both now
    # run one derivation, so a desalted container has to match a container
    # built from the desalted SMILES outright.
    built = MolContainer("CC(=O)C", "acetone", 0, {})
    desalted = MolContainer("CC(=O)C.CCO", "salt", 0, {})
    desalted.update_orig_smi("CC(=O)C")

    for field in (
        "orig_smi",
        "orig_smi_deslt",
        "orig_smi_canonical",
        "num_nonaro_rngs",
        "num_specif_chiral_cntrs",
        "num_unspecif_chiral_cntrs",
    ):
        assert getattr(desalted, field) == getattr(built, field), field


def test_update_orig_smi_stamps_the_container_index_on_the_rebuilt_mol() -> None:
    # The rebuilt reference molecule is handed to later steps (the ionization
    # fallback seeds a variant from it), and every step regroups its work by
    # contnr_idx, so the rebuild has to carry the container's index the way
    # construction does. update_orig_smi used to leave it at the MyMol default.
    contnr = MolContainer("CC(=O)C.CCO", "salt", 7, {})
    contnr.update_orig_smi("CC(=O)C")
    assert contnr.mol_orig_frm_inp_smi.contnr_idx == 7


def test_copy_of_orig_mol_is_independent_of_the_reference_mol() -> None:
    contnr = MolContainer("CCO", "ethanol", 2, {})
    contnr.mol_orig_frm_inp_smi.genealogy.append("CCO (source)")

    mol_copy = contnr.copy_of_orig_mol()

    assert mol_copy is not contnr.mol_orig_frm_inp_smi
    assert mol_copy.smiles() == contnr.mol_orig_frm_inp_smi.smiles()
    assert mol_copy.genealogy == ["CCO (source)"]
    assert mol_copy.contnr_idx == 2

    # The point of the copy: the working set can be rewritten in place without
    # touching the container's record of the input.
    mol_copy.genealogy.append("marker")
    mol_copy.make_first_3d_conf_no_min()

    assert contnr.mol_orig_frm_inp_smi.genealogy == ["CCO (source)"]
    assert contnr.mol_orig_frm_inp_smi.conformers == []


def test_update_idx_propagates_to_original_mol() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.update_idx(5)
    assert contnr.contnr_idx == 5
    assert contnr.mol_orig_frm_inp_smi.contnr_idx == 5


def test_update_idx_requires_int() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    with pytest.raises(Exception):
        contnr.update_idx("5")
