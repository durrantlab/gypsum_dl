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


def test_load_smiles_file_renames_duplicates(tmp_path) -> None:
    path = tmp_path / "input.smi"
    path.write_text("CCO\tethanol\nOCC\tethanol\nCCCO\tethanol\n")
    assert [d[1] for d in load_smiles_file(str(path))] == [
        "ethanol",
        "ethanol_copy_2",
        "ethanol_copy_3",
    ]


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


def test_load_sdf_file_skips_unparseable_records(tmp_path) -> None:
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", "ethanol")
    path = tmp_path / "input.sdf"
    _write_sdf(str(path), [mol])
    with open(path, "a") as f:
        f.write(BROKEN_SDF_RECORD)
    assert [d[1] for d in load_sdf_file(str(path))] == ["ethanol"]