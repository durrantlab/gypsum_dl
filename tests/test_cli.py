"""Tests for the command-line entry point, including the JSON override path."""

import json
import os
import sys

import pytest

from gypsum_dl import run


def test_cli_prepares_molecules_from_command_line(tmp_path, monkeypatch) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    output_folder = tmp_path / "cli_out"
    output_folder.mkdir()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "gypsum-dl",
            "--source",
            str(src),
            "--output_folder",
            str(output_folder),
            "--job_manager",
            "serial",
            "--max_variants_per_compound",
            "1",
            "--thoroughness",
            "1",
            "--2d_output_only",
        ],
    )
    run.main()
    assert os.path.exists(os.path.join(str(output_folder), "gypsum_dl_success.sdf"))


def test_cli_reads_parameters_from_json(
    tmp_path, monkeypatch, capsys: pytest.CaptureFixture[str]
) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    output_folder = tmp_path / "json_out"
    params_path = tmp_path / "params.json"
    params_path.write_text(
        json.dumps(
            {
                "source": str(src),
                "output_folder": str(output_folder),
                "job_manager": "serial",
                "max_variants_per_compound": 1,
                "thoroughness": 1,
                "2d_output_only": True,
            }
        )
    )
    monkeypatch.setattr(sys, "argv", ["gypsum-dl", "--json", str(params_path)])
    run.main()
    out = capsys.readouterr().out
    assert "overrides all other flags" in out
    assert os.path.exists(os.path.join(str(output_folder), "gypsum_dl_success.sdf"))


def test_cli_requires_a_source(monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["gypsum-dl", "--job_manager", "serial"])
    with pytest.raises(Exception, match="source"):
        run.main()