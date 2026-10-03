"""Unit tests for SMI and SDF input parsing, including the naming fallbacks."""

import pytest
from rdkit import Chem

from gypsum_dl.steps.io.LoadFiles import load_sdf_file, load_smiles_file

BROKEN_SDF_RECORD = """broken
     RDKit          2D

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0
M  END
$$$$
"""

# A record with a blank title line and no atoms. RDKit parses this into a valid
# (empty) Mol, so it reaches the naming code before anything notices it has no
# SMILES to contribute.
EMPTY_SDF_RECORD = """
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V2000
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
    assert props["activity"] == "1.5"


def test_load_sdf_file_keeps_property_values_verbatim(tmp_path) -> None:
    # Regression: properties were read with GetPropsAsDict, which converts
    # numeric-looking strings, so leading zeros, trailing zeros, and exponent
    # notation were lost before the values were written back to the output.
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", "ethanol")
    raw_values = {
        "CatalogID": "00123",
        "activity": "1.50",
        "conc": "3E4",
        "barcode": "12345678901234",
    }
    for key, val in raw_values.items():
        mol.SetProp(key, val)
    path = tmp_path / "input.sdf"
    _write_sdf(str(path), [mol])
    props = load_sdf_file(str(path))[0][2]
    assert props == raw_values


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


def test_load_sdf_file_skips_atomless_records_without_shifting_names(
    tmp_path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Regression: an atomless record was named, counted, and registered in
    # name_set before being dropped at the very end, so the untitled molecule
    # that follows it was named untitled_1_molnum_1 even though it is the first
    # molecule in the output. The record also vanished with no warning.
    untitled = Chem.MolFromSmiles("CCO")
    path = tmp_path / "input.sdf"
    path.write_text(EMPTY_SDF_RECORD + Chem.MolToMolBlock(untitled) + "$$$$\n")

    data = load_sdf_file(str(path))

    assert [d[0] for d in data] == ["CCO"]
    assert [d[1] for d in data] == ["untitled_0_molnum_0"]
    assert "Skipping an SDF record" in capsys.readouterr().out


def test_load_sdf_file_skips_atomless_records_without_claiming_names(
    tmp_path, capsys: pytest.CaptureFixture[str]
) -> None:
    # Same defect, seen from the duplicate-renaming side: the dropped record
    # claimed its title in name_set, so a later molecule with that title was
    # retitled against a name that appears nowhere in the output.
    lig = Chem.MolFromSmiles("CCO")
    lig.SetProp("_Name", "lig")
    path = tmp_path / "input.sdf"
    path.write_text(
        EMPTY_SDF_RECORD.replace("\n", "lig\n", 1) + Chem.MolToMolBlock(lig) + "$$$$\n"
    )

    assert [d[1] for d in load_sdf_file(str(path))] == ["lig"]
    assert "Skipping an SDF record" in capsys.readouterr().out


def test_load_sdf_file_skips_unparseable_records(tmp_path) -> None:
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("_Name", "ethanol")
    path = tmp_path / "input.sdf"
    _write_sdf(str(path), [mol])
    with open(path, "a") as f:
        f.write(BROKEN_SDF_RECORD)
    assert [d[1] for d in load_sdf_file(str(path))] == ["ethanol"]
