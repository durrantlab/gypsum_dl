"""Tests for the logging environment variables read at package import."""

import os
import subprocess
import sys
from pathlib import Path

import pytest


def _run_in_fresh_interpreter(
    code: str, **env_overrides: str
) -> subprocess.CompletedProcess[str]:
    """Run code in a fresh interpreter with the given environment.

    The logging variables are read once, when gypsum_dl is first imported, so
    their effect only shows up in a process that has not imported it yet.

    Args:
        code: Python source to run with `python -c`.
        **env_overrides: Environment variables to set for the child process.

    Returns:
        The completed subprocess, with stdout and stderr captured as text.
    """
    env = dict(os.environ)
    env.update(env_overrides)
    return subprocess.run(
        [sys.executable, "-c", code],
        env=env,
        capture_output=True,
        text=True,
    )


def _import_gypsum_dl(**env_overrides: str) -> subprocess.CompletedProcess[str]:
    """Import the package in a fresh interpreter with the given environment.

    Args:
        **env_overrides: Environment variables to set for the child process.

    Returns:
        The completed subprocess, with stdout and stderr captured as text.
    """
    return _run_in_fresh_interpreter("import gypsum_dl", **env_overrides)


@pytest.mark.parametrize("value", ["true", "yes", "on", "1", "True", "TRUE"])
def test_truthy_log_flags_do_not_break_the_import(value: str) -> None:
    # Regression: the flag was parsed with ast.literal_eval, which accepts only
    # Python literals, so every spelling but "True" raised ValueError from
    # inside ast while the package was still importing.
    result = _import_gypsum_dl(GYPSUM_DL_LOG=value)

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("value", ["false", "no", "off", "0", "False", ""])
def test_falsy_log_flags_do_not_break_the_import(value: str) -> None:
    result = _import_gypsum_dl(GYPSUM_DL_LOG=value)

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("level", ["DEBUG", "debug", "10", "20"])
def test_log_levels_accept_names_and_numbers(level: str) -> None:
    # Loguru takes either a level name or a numeric level, but the value was
    # forced through int(), so "DEBUG" crashed the import.
    result = _import_gypsum_dl(GYPSUM_DL_LOG="true", GYPSUM_DL_LOG_LEVEL=level)

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("value", ["true", "false", "True", "0"])
def test_stdout_flag_spellings_do_not_break_the_import(value: str) -> None:
    result = _import_gypsum_dl(GYPSUM_DL_LOG="true", GYPSUM_DL_STDOUT=value)

    assert result.returncode == 0, result.stderr


def test_unrecognized_flag_value_falls_back_to_the_default() -> None:
    # An unparseable value must not take the package down with it; logging
    # simply stays at its default setting.
    result = _import_gypsum_dl(GYPSUM_DL_LOG="maybe")

    assert result.returncode == 0, result.stderr


_LOG_MARKER_CODE = "from gypsum_dl import utils; utils.log('gypsum_log_marker')"


def test_log_file_receives_messages_without_color_codes(tmp_path: Path) -> None:
    """Messages must reach the log file, as plain text, and print only once.

    Regression: enable_logging configured loguru sinks, but every message went
    through print in utils.log, so GYPSUM_DL_LOG_FILE_PATH produced an empty
    file. The file sink was also colorized, which would have written ANSI
    escape codes into it.
    """
    log_path = tmp_path / "run.log"

    result = _run_in_fresh_interpreter(
        _LOG_MARKER_CODE,
        GYPSUM_DL_LOG="1",
        GYPSUM_DL_LOG_FILE_PATH=str(log_path),
    )

    assert result.returncode == 0, result.stderr
    contents = log_path.read_text()
    assert "gypsum_log_marker" in contents
    assert "\x1b[" not in contents
    # The stdout sink is off by default, so the print is the only copy.
    assert result.stdout.count("gypsum_log_marker") == 1


def test_messages_stay_out_of_loguru_when_logging_is_off() -> None:
    """With logging off, messages are printed and loguru stays silent.

    loguru's default handler writes to stderr, so a message that escaped
    logger.disable would show up there.
    """
    result = _run_in_fresh_interpreter(_LOG_MARKER_CODE, GYPSUM_DL_LOG="0")

    assert result.returncode == 0, result.stderr
    assert result.stdout.count("gypsum_log_marker") == 1
    assert "gypsum_log_marker" not in result.stderr
