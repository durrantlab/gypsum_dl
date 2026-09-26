"""Unit tests for MyMol and MyConformer."""

import pytest
from rdkit import Chem

from gypsum_dl import MyMol


def test_mymol_from_rdkit_mol_sets_canonical_smiles() -> None:
    mol = MyMol.MyMol(Chem.MolFromSmiles("OCC"), "ethanol")
    assert mol.can_smi == "CCO"
    assert mol.name == "ethanol"


def test_mymol_from_rdkit_mol_survives_smiles_conversion_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: when MolToSmiles raises for an RDKit-mol starter, __init__
    # used to reference the unbound local `smiles` and die with
    # UnboundLocalError instead of falling back to can_smi=False.
    def boom(*args, **kwargs):
        raise ValueError("cannot canonicalize")

    monkeypatch.setattr(MyMol.Chem, "MolToSmiles", boom)
    mol = MyMol.MyMol(Chem.MolFromSmiles("OCC"), "ethanol")
    assert mol.can_smi is False
    assert mol.orig_smi == ""


def test_smiles_noh_strips_explicit_hydrogens() -> None:
    assert MyMol.MyMol("CCO").smiles(True) == "CCO"


def test_smiles_is_cached() -> None:
    mol = MyMol.MyMol("CCO")
    assert mol.smiles() == "CCO"
    assert mol.smiles() == "CCO"
    assert mol.smiles(True) == mol.smiles(True)


def test_comparison_operators_follow_canonical_smiles() -> None:
    a = MyMol.MyMol("CCO")
    b = MyMol.MyMol("OCC")
    c = MyMol.MyMol("CCC")
    assert a == b
    assert a != c
    assert a <= b
    assert a >= b
    assert (a < c) == (hash(a) < hash(c))
    assert (a > c) == (hash(a) > hash(c))
    assert a is not None


def test_standardize_smiles_is_cached() -> None:
    mol = MyMol.MyMol("CCO")
    first = mol.standardize_smiles()
    assert first == mol.standardize_smiles()


def test_count_hyd_bnd_to_carb() -> None:
    assert MyMol.MyMol("CCO").count_hyd_bnd_to_carb() == 5


def test_get_idxs_of_nonaro_rng_atms_is_cached() -> None:
    mol = MyMol.MyMol("C1CCCCC1")
    rings = mol.get_idxs_of_nonaro_rng_atms()
    assert len(rings) == 1
    assert mol.get_idxs_of_nonaro_rng_atms() is rings


def test_get_idxs_of_nonaro_rng_atms_ignores_aromatic_rings() -> None:
    assert MyMol.MyMol("c1ccccc1").get_idxs_of_nonaro_rng_atms() == []


def test_chiral_center_helpers_are_cached() -> None:
    mol = MyMol.MyMol("CC(N)C(=O)O")
    unassigned = mol.chiral_cntrs_w_unasignd()
    assert len(unassigned) == 1
    assert mol.chiral_cntrs_w_unasignd() is unassigned
    assigned = mol.chiral_cntrs_only_asignd()
    assert assigned == []
    assert mol.chiral_cntrs_only_asignd() is assigned


def test_get_double_bonds_without_stereochemistry_finds_unspecified_bond() -> None:
    assert len(MyMol.MyMol("CC=CC").get_double_bonds_without_stereochemistry()) == 1


def test_get_double_bonds_without_stereochemistry_ignores_specified_bond() -> None:
    mol = MyMol.MyMol(Chem.MolFromSmiles(r"C/C=C/C"))
    assert mol.get_double_bonds_without_stereochemistry() == []


def test_remove_bizarre_substruc_flags_carbanion() -> None:
    mol = MyMol.MyMol("CC[CH2-]")
    assert mol.remove_bizarre_substruc() is True
    assert mol.remove_bizarre_substruc() is True


def test_remove_bizarre_substruc_allows_normal_molecule() -> None:
    mol = MyMol.MyMol("CCO")
    assert mol.remove_bizarre_substruc() is False
    assert mol.remove_bizarre_substruc() is False


def test_remove_bizarre_substruc_survives_non_string_can_smi() -> None:
    # Regression (M8): can_smi is False after a failed MolToSmiles and None
    # after a failed smiles(); `s in self.can_smi` then raised TypeError
    # ("argument of type 'bool'/'NoneType' is not iterable"), which became a
    # hang under multiprocessing. The method must return a bool instead.
    mol = MyMol.MyMol("CCO")
    mol.can_smi = False
    result = mol.remove_bizarre_substruc()
    assert isinstance(result, bool)

    mol2 = MyMol.MyMol("CCO")
    mol2.can_smi = None
    assert isinstance(mol2.remove_bizarre_substruc(), bool)


def test_get_frags_of_orig_smi_single_fragment_returns_self() -> None:
    mol = MyMol.MyMol("CCO")
    assert mol.get_frags_of_orig_smi() == [mol]


def test_get_frags_of_orig_smi_splits_salts() -> None:
    assert len(MyMol.MyMol("CCO.CC").get_frags_of_orig_smi()) == 2


def test_inherit_contnr_props() -> None:
    source = MyMol.MyMol("CCO", "ethanol")
    source.contnr_idx = 3
    target = MyMol.MyMol("CCC")
    target.inherit_contnr_props(source)
    assert target.contnr_idx == 3
    assert target.name == "ethanol"
    assert target.orig_smi == source.orig_smi


def test_set_all_rdkit_mol_props_records_genealogy() -> None:
    mol = MyMol.MyMol("CCO", "ethanol")
    mol.mol_props["activity"] = 1.5
    mol.genealogy = ["CCO (source)"]
    mol.set_all_rdkit_mol_props()
    assert mol.rdkit_mol.GetProp("SMILES") == "CCO"
    assert mol.rdkit_mol.GetProp("activity") == "1.5"
    assert mol.rdkit_mol.GetProp("Genealogy") == "CCO (source)"
    assert mol.rdkit_mol.GetProp("_Name") == "ethanol"


def test_make_first_3d_conf_no_min_is_idempotent() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    assert len(mol.conformers) == 1
    mol.make_first_3d_conf_no_min()
    assert len(mol.conformers) == 1


def test_make_first_3d_conf_no_min_survives_failed_reprotanation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression (M13): try_reprotanation can return None. That None was
    # assigned straight to rdkit_mol, and add_conformers -> MyConformer then did
    # copy.deepcopy(None).RemoveAllConformers(), raising AttributeError inside
    # pick_lowest_enrgy_mols. The method must bail out quietly instead.
    monkeypatch.setattr(MyMol.MOH, "try_reprotanation", lambda mol: None)
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    assert mol.conformers == []


def test_myconformer_with_none_rdkit_mol_is_marked_failed() -> None:
    # Regression (M13): constructing a MyConformer from a MyMol whose rdkit_mol
    # is None must not crash (deepcopy(None).RemoveAllConformers()). It should
    # instead flag itself failed (mol is False) so add_conformers skips it.
    mol = MyMol.MyMol("CCO")
    mol.rdkit_mol = None
    conf = MyMol.MyConformer(mol)
    assert conf.mol is False


def test_add_conformers_sorts_by_energy() -> None:
    mol = MyMol.MyMol("CCCCCC")
    # `MyConformer.rmsd_to_me` rebuilds the molecule from SMILES and
    # reprotonates it, so it only matches conformers whose parent already
    # carries explicit hydrogens. The pipeline guarantees this by calling
    # `make_first_3d_conf_no_min` before any RMSD-based pruning; do the same
    # here rather than embedding an implicit-hydrogen molecule.
    mol.make_first_3d_conf_no_min()
    mol.add_conformers(3, 0.1, True)
    energies = [conf.energy for conf in mol.conformers]
    assert len(energies) >= 1
    assert energies == sorted(energies)


def test_load_conformers_into_rdkit_mol() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    mol.load_conformers_into_rdkit_mol()
    assert mol.rdkit_mol.GetNumConformers() == 1


def test_conformer_coords_and_energy() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    assert conf.coords().shape[0] == conf.mol.GetNumAtoms()
    assert isinstance(conf.get_energy(), float)


def test_conformer_minimize_is_idempotent() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    conf.minimize()
    energy = conf.energy
    conf.minimize()
    assert conf.energy == energy


def test_conformer_write_pdb_file(tmp_path) -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    out = tmp_path / "conf.pdb"
    mol.conformers[0].write_pdb_file(str(out))
    # Small molecules carry no residue information, so RDKit emits HETATM.
    assert "HETATM" in out.read_text()


def test_conformer_rmsd_between_identical_conformers_is_zero() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    duplicate = MyMol.MyConformer(mol, conf.conformer())
    assert conf.rmsd_to_me(duplicate) == pytest.approx(0.0, abs=1e-6)


def test_coord_3d_err_warning_is_logged(capsys: pytest.CaptureFixture[str]) -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    mol.conformers[0].coord_3d_err_warning(None)
    assert "WARNING" in capsys.readouterr().out