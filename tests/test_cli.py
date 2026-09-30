"""Tests for the command-line entry point, including the JSON override path."""

import json
import os
import re
import runpy
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
    # (The warning only makes sense when some other flag is actually supplied
    # on the command line; --num_processors no longer leaks in via a non-None
    # argparse default.)
    monkeypatch.setattr(
        sys,
        "argv",
        ["gypsum-dl", "--json", str(params_path), "--num_processors", "1"],
    )
    run.main()
    out = capsys.readouterr().out
    assert "overrides all other flags" in out
    assert os.path.exists(os.path.join(str(output_folder), "gypsum_dl_success.sdf"))


@pytest.mark.parametrize(
    "extra_args",
    [
        ["--use_durrant_lab_filters"],
        ["--job_manager", "mpi"],
    ],
    ids=["store_true_flag", "job_manager"],
)
def test_cli_warns_when_json_overrides_any_argument(
    tmp_path,
    monkeypatch,
    capsys: pytest.CaptureFixture[str],
    extra_args: list[str],
) -> None:
    # Regression (F1): the override warning fired only for eight named
    # parameters, and --job_manager plus the fourteen store_true flags could
    # never fire it at all because argparse always materialized them. So
    # `--json p.json --job_manager mpi` ran one multiprocessing job per rank,
    # all writing to the same output folder, with nothing said about it.
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
    monkeypatch.setattr(
        sys, "argv", ["gypsum-dl", "--json", str(params_path)] + extra_args
    )
    run.main()
    out = capsys.readouterr().out
    assert "overrides all other flags" in out
    # The json file still wins: the run has to be the one it describes.
    assert os.path.exists(os.path.join(str(output_folder), "gypsum_dl_success.sdf"))


def test_cli_records_only_the_arguments_the_user_supplied(monkeypatch) -> None:
    # Regression (F1): --job_manager carried a non-None argparse default and
    # the store_true flags defaulted to False, so both survived the
    # None-stripping in main() and reached prepare_molecules as though the
    # user had typed them. That is what made the discarded-argument warning
    # impossible to compute, and it also let argparse, rather than
    # set_parameters, own the defaults.
    captured: dict = {}
    monkeypatch.setattr(run, "prepare_molecules", lambda args: captured.update(args))
    monkeypatch.setattr(sys, "argv", ["gypsum-dl", "--source", "x.smi"])
    run.main()
    assert captured == {"source": "x.smi"}


def test_cli_passes_a_supplied_store_true_flag_through(monkeypatch) -> None:
    captured: dict = {}
    monkeypatch.setattr(run, "prepare_molecules", lambda args: captured.update(args))
    monkeypatch.setattr(
        sys,
        "argv",
        ["gypsum-dl", "--source", "x.smi", "--use_durrant_lab_filters"],
    )
    run.main()
    assert captured == {"source": "x.smi", "use_durrant_lab_filters": True}


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


def test_cli_rejects_the_removed_cache_prerun_flag(monkeypatch) -> None:
    # --cache_prerun only ever short-circuited the run so that one rank wrote
    # the __pycache__ files before the rest started. Nothing in this fork read
    # it, so the documented `mpirun -n 1 ... -c` line reached prepare_molecules
    # with no source and died there. The warm-up now lives in the MPI job
    # script itself (see tests/files/Pitt_CRC/mpi_gypsum_dl.sh), so the flag is
    # gone rather than reimplemented.
    monkeypatch.setattr(sys, "argv", ["gypsum-dl", "-c"])
    with pytest.raises(SystemExit):
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


README = Path(__file__).resolve().parent.parent / "README.md"


def _parse_option_help(block: str) -> dict[str, str]:
    """Map each option to its help text, as rendered by argparse.

    Both the README block and the real --help output are argparse renderings of
    the same parser, so parsing them the same way lets them be compared without
    depending on the terminal width (which only changes where lines wrap) or on
    the argparse version (which only changes how invocations are joined). The
    invocation is separated from its help by at least two spaces, and never
    contains two consecutive spaces itself, so that gap is what splits them.

    Args:
        block: The body of an argparse options section.

    Returns:
        Help text, whitespace-collapsed, keyed by the option's first long form.
    """
    entries: dict[str, str] = {}
    key = ""
    for line in block.split("\n"):
        if not line.strip():
            continue
        if line.startswith("  ") and line[2:3] == "-":
            parts = re.split(r"\s{2,}", line.strip(), maxsplit=1)
            invocation = parts[0]
            long_forms = [
                token
                for token in re.split(r"[\s,=]+", invocation)
                if token.startswith("--")
            ]
            key = long_forms[0] if long_forms else invocation
            entries[key] = parts[1] if len(parts) > 1 else ""
        elif key:
            entries[key] += " " + line.strip()
    return {option: " ".join(text.split()) for option, text in entries.items()}


def _readme_option_block() -> str:
    """Extract the fenced command-line parameter block from the README.

    Returns:
        The text inside the first ```text fence, which documents the options.
    """
    text = README.read_text(encoding="utf-8")
    blocks = re.findall(r"^```text\n(.*?)^```", text, flags=re.MULTILINE | re.DOTALL)
    assert blocks, "README.md has no ```text command-line parameter block."
    return blocks[0]


def _cli_option_block(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> str:
    """Run the CLI with --help and return just its options section.

    Args:
        monkeypatch: Fixture used to fix the terminal width and argv.
        capsys: Fixture used to capture the help output.

    Returns:
        The options section of the help output, up to the epilog.
    """
    monkeypatch.setenv("COLUMNS", "80")
    monkeypatch.setattr(sys, "argv", ["gypsum-dl", "--help"])
    with pytest.raises(SystemExit):
        run.main()
    out = capsys.readouterr().out
    sections = re.split(
        r"^(?:options|optional arguments):\n", out, flags=re.MULTILINE, maxsplit=1
    )
    assert len(sections) == 2, "argparse help output has no options section."
    return sections[1].split("\n\n")[0]


def test_readme_documents_the_help_text_the_cli_prints(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    # Regression: the README credited --add_html_output with opening a browser,
    # which nothing does, and omitted the --num_processors default. The README
    # block is a copy of --help, so it can only be trusted if it is pinned.
    assert _parse_option_help(_readme_option_block()) == _parse_option_help(
        _cli_option_block(monkeypatch, capsys)
    )


def test_package_runs_as_a_module(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    # Regression: the docs launched mpi mode through run_gypsum_dl.py, which a
    # pip install does not provide, and the package had no __main__, so
    # "python -m mpi4py -m gypsum_dl" had nothing to run.
    monkeypatch.setattr(sys, "argv", ["gypsum_dl", "--help"])
    with pytest.raises(SystemExit):
        runpy.run_module("gypsum_dl", run_name="__main__")
    assert "EXAMPLES OF USE" in capsys.readouterr().out
