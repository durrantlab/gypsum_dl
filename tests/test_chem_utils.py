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
        mol.make_first_3d_conf_no_min = lambda m=mol: None  # no-op; conformers already set

    # thoroughness=3 so random_sample draws all 3 candidates; run 20 times to
    # defeat the internal shuffle.
    for _ in range(20):
        kept = chem_utils.pick_lowest_enrgy_mols(mols, 1, 3)
        assert len(kept) == 1
        assert kept[0].smiles() == "CCC", (
            f"Expected lowest-energy mol 'CCC' but got '{kept[0].smiles()}'"
        )


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
