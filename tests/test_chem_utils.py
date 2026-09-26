"""Unit tests for the variant-selection helpers in chem_utils."""

from gypsum_dl import MyMol, chem_utils
from gypsum_dl.MolContainer import MolContainer


def test_uniq_mols_in_list_removes_duplicates() -> None:
    mols = [MyMol.MyMol("CCO"), MyMol.MyMol("OCC"), MyMol.MyMol("CCC")]
    assert len(chem_utils.uniq_mols_in_list(mols)) == 2


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