"""Tests for the command-line entry point, including the JSON override path."""

import json
import os
import sys
from pathlib import Path

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
    # Pass an overridable flag alongside --json so the override warning fires.
    # (The warning only makes sense when a json_warning_list flag is actually
    # supplied on the command line; --num_processors no longer leaks in via a
    # non-None argparse default.)
    monkeypatch.setattr(
        sys,
        "argv",
        ["gypsum-dl", "--json", str(params_path), "--num_processors", "1"],
    )
    run.main()
    out = capsys.readouterr().out
    assert "overrides all other flags" in out
    assert os.path.exists(os.path.join(str(output_folder), "gypsum_dl_success.sdf"))


def test_cli_num_processors_defaults_to_all_cores(monkeypatch) -> None:
    # Regression (B18): the argparse default for --num_processors was 1, which
    # (unlike the other numeric flags) is never None and so never stripped,
    # pinning the CLI to a single core regardless of set_parameters' -1 default.
    # With the default now None, omitting the flag lets set_parameters decide.
    captured: dict = {}
    monkeypatch.setattr(run, "prepare_molecules", lambda args: captured.update(args))
    monkeypatch.setattr(sys, "argv", ["gypsum-dl", "--source", "x.smi"])
    run.main()
    assert "num_processors" not in captured


def test_cli_requires_a_source(monkeypatch) -> None:
    monkeypatch.setattr(sys, "argv", ["gypsum-dl", "--job_manager", "serial"])
    with pytest.raises(Exception, match="source"):
        run.main()


def test_cli_accepts_random_seed(monkeypatch) -> None:
    # Regression: set_parameters and seed_random_number_generators both
    # understood random_seed, but argparse did not, so the only way to reach it
    # was a JSON parameter file, which overrides every other flag. Anyone
    # wanting a reproducible run plus command-line control could not have both.
    captured: dict = {}
    monkeypatch.setattr(run, "prepare_molecules", lambda args: captured.update(args))
    monkeypatch.setattr(
        sys, "argv", ["gypsum-dl", "--source", "x.smi", "--random_seed", "42"]
    )
    run.main()
    assert captured["random_seed"] == 42


def test_cli_omits_random_seed_when_not_given(monkeypatch) -> None:
    # The flag must not carry a non-None argparse default, or it would survive
    # the None-stripping in main() and shadow set_parameters' own default.
    captured: dict = {}
    monkeypatch.setattr(run, "prepare_molecules", lambda args: captured.update(args))
    monkeypatch.setattr(sys, "argv", ["gypsum-dl", "--source", "x.smi"])
    run.main()
    assert "random_seed" not in captured


def _sdf_molecule_records(path: Path) -> list[str]:
    """Read an SDF's molecule records, dropping the leading parameter record.

    save_to_sdf writes the run's parameters as the first record, and those
    carry a wall-clock start_time, so comparing two runs' whole files could
    never show them as identical however well seeded the runs were.

    Args:
        path: Path to the SDF file to read.

    Returns:
        The records describing actual molecules, in file order.
    """
    return path.read_text().split("$$$$")[1:]


def test_cli_random_seed_makes_a_serial_run_reproducible(tmp_path, monkeypatch) -> None:
    # The flag has to actually reach the generators, not merely parse: this is
    # the behavior the parameter exists for.
    src = tmp_path / "input.smi"
    src.write_text("CC(N)C(=O)O\talanine\n")

    records = []
    for run_name in ("first", "second"):
        output_folder = tmp_path / run_name
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
                "2",
                "--thoroughness",
                "1",
                "--random_seed",
                "202609",
                "--2d_output_only",
            ],
        )
        run.main()
        records.append(_sdf_molecule_records(output_folder / "gypsum_dl_success.sdf"))

    assert records[0] == records[1]
