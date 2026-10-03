"""Simplify your molecular simulation workflow."""

from typing import Any

import os
import sys

from loguru import logger

logger.disable("gypsum_dl")

LOG_FORMAT = (
    "<green>{time:HH:mm:ss}</green> | "
    "<level>{level: <8}</level> | "
    "<cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>"
)


def enable_logging(
    level_set: int | str,
    stdout_set: bool = False,
    file_path: str | None = None,
    log_format: str = LOG_FORMAT,
) -> None:
    r"""Enable logging.

    Args:
        level_set: Requested log level: `10` is debug, `20` is info.
        stdout_set: Also send log records to stdout. Off by default because
            every message is already printed there.
        file_path: Also write logs to files here.
        log_format: The loguru format string for each record.
    """
    config: dict[str, Any] = {"handlers": []}
    if stdout_set:
        config["handlers"].append(
            {
                "sink": sys.stdout,
                "level": level_set,
                "format": log_format,
                "colorize": True,
            }
        )
    if isinstance(file_path, str):
        config["handlers"].append(
            {
                "sink": file_path,
                "level": level_set,
                "format": log_format,
                "colorize": False,
            }
        )
    # https://loguru.readthedocs.io/en/stable/api/logger.html#loguru._logger.Logger.configure
    logger.configure(**config)

    logger.enable("gypsum_dl")


_TRUE_STRINGS = frozenset({"1", "true", "yes", "on"})
_FALSE_STRINGS = frozenset({"0", "false", "no", "off", ""})


def _env_flag(name: str, default: bool) -> bool:
    """Read a boolean environment variable without failing the package import.

    `ast.literal_eval` accepts only Python literals, so the ordinary spellings a
    user reaches for raised `ValueError` while `gypsum_dl` was still importing.

    Args:
        name: Environment variable to read.
        default: Value to use when the variable is unset or unrecognized.

    Returns:
        The parsed flag.
    """
    raw = os.environ.get(name)
    if raw is None:
        return default
    normalized = raw.strip().lower()
    if normalized in _TRUE_STRINGS:
        return True
    if normalized in _FALSE_STRINGS:
        return False
    return default


def _env_log_level(default: int) -> int | str:
    """Read a loguru level from the environment, accepting names or numbers.

    Loguru takes either form, so rejecting `DEBUG` with a bare `int()` call
    turned a reasonable setting into an import-time crash.

    Args:
        default: Numeric level to use when the variable is unset.

    Returns:
        A numeric level, or an uppercased loguru level name.
    """
    raw = os.environ.get("GYPSUM_DL_LOG_LEVEL")
    if raw is None:
        return default
    raw = raw.strip()
    try:
        return int(raw)
    except ValueError:
        return raw.upper()


if _env_flag("GYPSUM_DL_LOG", False):
    enable_logging(
        _env_log_level(20),
        _env_flag("GYPSUM_DL_STDOUT", False),
        os.environ.get("GYPSUM_DL_LOG_FILE_PATH", None),
    )
