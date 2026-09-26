"""Unit tests for SMI and SDF input parsing, including the naming fallbacks."""

from rdkit import Chem

from gypsum_dl.steps.io.LoadFiles import load_sdf_file, load_smiles_file

BROKEN_SDF_RECORD = """broken
     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
M  END
$$$$
"""


def _write_sdf(path: str, mols: list) -> None:
    """Write RDKit molecules to an SDF file.

    Args:
        path: Destination file path.
        mols: RDKit molecules to serialize.
    """
    writer = Chem.SDWriter(path)
    for mol in mols:
        writer.write(mol)
    writer.close()


def test_load_smiles_file_reads_names_and_skips_blank_lines(tmp_path) -> None:
    path = tmp_path / "input.smi"
    path.write_text("CCO\tethanol\n\nCCC propane\n")
    data = load_smiles_file(str(path))
    assert [d[0] for d in data] == ["CCO", "CCC"]
    assert [d[1] for d in data] == ["ethanol", "propane"]


def test_load_smiles_file_names_untitled_entries(tmp_path) -> None:
    path = tmp_path / "input.smi"
    path.write_text("CCO\nCCC\n")
    assert [d[1] for d in load_smiles_file(str(path))] == [
        "untitled_line_1",
        "untitled_line_2",
    ]


def test_load_smiles_file_line_numbers_account_for_blank_lines(tmp_path) -> None:
    # Regression (B15): the untitled-ligand name (and log line numbers) tracked
    # a counter that only advanced on non-blank lines, so blank lines threw the
    # reported line number off. The molecule below sits on file line 3.
    path = tmp_path / "input.smi"
    path.write_text("\n\nCCO\n")
    assert [d[1] for d in load_smiles_file(str(path))] == ["untitled_line_3"]


def test_load_smiles_file_reads_utf8_names(tmp_path) -> None:
    # Regression (B16): the file was opened with the platform default encoding,
    # which raises UnicodeDecodeError on a UTF-8 ligand name under a non-UTF-8
    # locale. It must be read as UTF-8.
    path = tmp_path / "input.smi"
    path.write_text("CCO café\n", encoding="utf-8")
    assert [d[1] for d in load_smiles_file(str(path))] == ["café"]


def test_load_smiles_file_renames_duplicates(tmp_path) -> None:
    path = tmp_path / "input.smi"
    path.write_text("CCO\tethanol\nOCC\tethanol\nCCCO\tethanol\n")
    assert [d[1] for d in load_smiles_file(str(path))] == [
        "ethanol",
        "ethanol_copy_2",
        "ethanol_copy_3",
    ]


def test_load_smiles_file_duplicate_renaming_avoids_existing_names(tmp_path) -> None:
    # Regression (B13): the generated "_copy_N" name was never checked against
    # the names already seen, so an input that itself contains "lig_copy_2"
    # produced two records with that same name.
    path = tmp_path / "input.smi"
    path.write_text("CCO lig_copy_2\nCCC lig\nCCCC lig\n")
    names = [d[1] for d in load_smiles_file(str(path))]
    assert names == ["lig_copy_2", "lig", "lig_copy_3"]
    assert len(set(names)) == 3


def test_load_sdf_file_reads_names_and_properties(tmp_path) -> None:
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", "ethanol")
    mol.SetProp("activity", "1.5")
    path = tmp_path / "input.sdf"
    _write_sdf(str(path), [mol])
    data = load_sdf_file(str(path))
    assert len(data) == 1
    smiles, name, props = data[0]
    assert smiles == "CCO"
    assert name == "ethanol"
    assert props["activity"] == 1.5


def test_load_sdf_file_names_untitled_molecules(tmp_path) -> None:
    path = tmp_path / "input.sdf"
    _write_sdf(str(path), [Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("CCC")])
    assert [d[1] for d in load_sdf_file(str(path))] == [
        "untitled_0_molnum_0",
        "untitled_1_molnum_1",
    ]


def test_load_sdf_file_renames_named_duplicates(tmp_path) -> None:
    # Regression: named duplicates must be deduplicated and mol_obj_counter must
    # advance for every record (not just untitled ones).
    lig1 = Chem.MolFromSmiles("CCO")
    lig1.SetProp("_Name", "lig")
    lig2 = Chem.MolFromSmiles("OCC")
    lig2.SetProp("_Name", "lig")
    untitled = Chem.MolFromSmiles("CCC")
    path = tmp_path / "input.sdf"
    _write_sdf(str(path), [lig1, lig2, untitled])
    assert [d[1] for d in load_sdf_file(str(path))] == [
        "lig",
        "lig_copy_2",
        "untitled_0_molnum_2",
    ]


def test_load_sdf_file_duplicate_renaming_avoids_existing_names(tmp_path) -> None:
    # Regression (B13): same collision as the SMI case above.
    mols = []
    for smiles, name in (("CCO", "lig_copy_2"), ("CCC", "lig"), ("CCCC", "lig")):
        mol = Chem.MolFromSmiles(smiles)
        mol.SetProp("_Name", name)
        mols.append(mol)
    path = tmp_path / "input.sdf"
    _write_sdf(str(path), mols)
    names = [d[1] for d in load_sdf_file(str(path))]
    assert names == ["lig_copy_2", "lig", "lig_copy_3"]
    assert len(set(names)) == 3


def test_load_sdf_file_skips_unparseable_records(tmp_path) -> None:
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", "ethanol")
    path = tmp_path / "input.sdf"
    _write_sdf(str(path), [mol])
    with open(path, "a") as f:
        f.write(BROKEN_SDF_RECORD)
    assert [d[1] for d in load_sdf_file(str(path))] == ["ethanol"]
