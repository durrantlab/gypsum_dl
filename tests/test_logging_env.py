"""Tests for the logging environment variables read at package import."""

import os
import subprocess
import sys

import pytest


def _import_gypsum_dl(**env_overrides: str) -> subprocess.CompletedProcess[str]:
    """Import the package in a fresh interpreter with the given environment.

    The variables in question are read once, at module scope, so the failure
    only shows up in a process that has not imported gypsum_dl yet.

    Args:
        **env_overrides: Environment variables to set for the child process.

    Returns:
        The completed subprocess, with stdout and stderr captured as text.
    """
    env = dict(os.environ)
    env.update(env_overrides)
    return subprocess.run(
        [sys.executable, "-c", "import gypsum_dl"],
        env=env,
        capture_output=True,
        text=True,
    )


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
