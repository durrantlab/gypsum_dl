"""Unit tests for parameter validation and the small helpers in start.py."""

import json
import os

import pytest

from gypsum_dl import start
from gypsum_dl.MolContainer import MolContainer


def test_detect_unassigned_bonds_accepts_valid_smiles() -> None:
    assert start.detect_unassigned_bonds("CCO") == "CCO"


def test_detect_unassigned_bonds_rejects_garbage() -> None:
    assert start.detect_unassigned_bonds("moosedogfacecat") is None


def test_make_type_dict_maps_scalar_types() -> None:
    type_dict = start.make_type_dict({"a": 1, "b": 1.5, "c": True, "d": "x"})
    assert type_dict == {"a": int, "b": float, "c": bool, "d": str}


def test_make_type_dict_rejects_unsupported_types() -> None:
    with pytest.raises(Exception, match="No value can have type"):
        start.make_type_dict({"a": None})


def test_merge_parameters_promotes_int_to_float() -> None:
    default = {"min_ph": 6.4}
    start.merge_parameters(default, {"min_ph": 7})
    assert isinstance(default["min_ph"], float)
    assert default["min_ph"] == pytest.approx(7.0)


def test_merge_parameters_rejects_unknown_parameter() -> None:
    with pytest.raises(Exception, match="Unrecognized parameter"):
        start.merge_parameters({"min_ph": 6.4}, {"bogus": 1})


def test_merge_parameters_rejects_wrong_type() -> None:
    with pytest.raises(Exception, match="must be of"):
        start.merge_parameters({"min_ph": 6.4}, {"min_ph": "high"})


def test_set_parameters_lowercases_keys_and_fills_defaults(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params = start.set_parameters({"SOURCE": str(src), "job_manager": "SERIAL"})
    assert params["source"] == os.path.abspath(str(src))
    assert params["job_manager"] == "serial"
    assert params["thoroughness"] == 3
    assert params["max_variants_per_compound"] == 5
    # Regression (M10): a random_seed parameter must exist so runs can be made
    # reproducible; the default (-1) leaves the RNG unseeded.
    assert params["random_seed"] == -1


def test_set_parameters_accepts_random_seed(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params = start.set_parameters({"source": str(src), "random_seed": 42})
    assert params["random_seed"] == 42


def test_finalize_params_requires_source() -> None:
    with pytest.raises(Exception, match="source"):
        start.finalize_params({"source": "", "output_folder": "./"})


def test_finalize_params_defaults_output_folder_next_to_source(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params = start.finalize_params(
        {
            "source": str(src),
            "output_folder": "",
            "add_pdb_output": False,
            "separate_output_files": False,
            "job_manager": "Serial",
        }
    )
    assert params["output_folder"].endswith(f"output{os.sep}")
    assert params["job_manager"] == "serial"


def test_finalize_params_derives_source_dir_by_dirname(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression (M11): source_dir used str.strip(basename), which strips a
    # character *set* off both ends rather than removing a suffix. On POSIX
    # absolute paths the leading "/" and the last "/" hide the bug, but any
    # path whose leading directory character also appears in the basename gets
    # mangled (e.g. "smiles_dir/mol.smi" -> "es_dir/", and on Windows the drive
    # letter is eaten). Pin abspath to identity so a relative path reaches the
    # source_dir computation intact, then assert it equals os.path.dirname.
    monkeypatch.setattr(start.os.path, "abspath", lambda p: p)
    params = start.finalize_params(
        {
            "source": "smiles_dir/mol.smi",
            "output_folder": "",
            "add_pdb_output": False,
            "separate_output_files": False,
            "job_manager": "serial",
        }
    )
    assert params["output_folder"] == "smiles_dir" + os.sep + "output" + os.sep


def test_add_mol_id_props_assigns_unique_ids() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    start.add_mol_id_props([contnr])
    assert contnr.mols[0].rdkit_mol.GetProp("UniqueID") == "1"


def test_deal_with_failed_molecules_writes_failure_file(tmp_path) -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    start.deal_with_failed_molecules([contnr], {"output_folder": str(tmp_path)})
    assert "ethanol" in (tmp_path / "gypsum_dl_failed.smi").read_text()


def test_deal_with_failed_molecules_skips_when_nothing_failed(tmp_path) -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    start.deal_with_failed_molecules([contnr], {"output_folder": str(tmp_path)})
    assert not os.path.exists(os.path.join(str(tmp_path), "gypsum_dl_failed.smi"))


def test_prepare_molecules_rejects_malformed_json(tmp_path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    with pytest.raises(Exception, match="properly formed"):
        start.prepare_molecules({"json": str(bad)})


def test_prepare_molecules_rejects_unknown_json_parameter(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps({"source": str(src), "bogus_flag": True}))
    with pytest.raises(Exception, match="Unrecognized parameter"):
        start.prepare_molecules({"json": str(params_path)})