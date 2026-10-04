"""Refuse to publish unless dist/ holds exactly the release that HEAD is tagged as.

setuptools-scm takes the version from git, so building from an untagged or
modified checkout quietly produces something like 2.0.1.dev3, and twine will
upload whatever is in dist/. This runs before every upload so that a missing
tag, a dirty tree, or leftover files from an earlier build stop the upload
instead of reaching PyPI.

Usage:
    python scripts/check_release.py
"""

import subprocess
import sys
from pathlib import Path
from typing import NoReturn

DIST = Path("dist")


def fail(message: str) -> NoReturn:
    """Print a reason and exit nonzero, so the pixi task chain stops.

    Args:
        message: What is wrong with the release.
    """

    sys.exit(f"Release check failed: {message}")


def dist_version() -> str:
    """Return the one version that every file in dist/ carries.

    Returns:
        The version shared by the single wheel and the single sdist.
    """

    files = sorted(p.name for p in DIST.glob("*")) if DIST.is_dir() else []
    wheels = [f for f in files if f.endswith(".whl")]
    sdists = [f for f in files if f.endswith(".tar.gz")]
    if len(wheels) != 1 or len(sdists) != 1 or len(files) != 2:
        fail(f"expected one wheel and one sdist in dist/, found {files}")

    # Wheel names are name-version-tags.whl; sdist names are name-version.tar.gz.
    wheel_version = wheels[0].split("-")[1]
    sdist_version = sdists[0][: -len(".tar.gz")].rsplit("-", 1)[1]
    if wheel_version != sdist_version:
        fail(f"wheel is {wheel_version} but sdist is {sdist_version}")
    return wheel_version


def head_tag() -> str:
    """Return the tag on the current commit.

    Returns:
        The tag name, such as v2.0.0.
    """

    result = subprocess.run(
        ["git", "describe", "--tags", "--exact-match", "HEAD"],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        fail("the current commit has no tag; run git tag vX.Y.Z first")
    return result.stdout.strip()


def main() -> None:
    """Check the built version against the tag and report what will upload."""

    version = dist_version()
    if "dev" in version or "+" in version:
        fail(
            f"built version is {version}, not a release; commit or stash any "
            "changes, make sure HEAD is tagged, and rebuild"
        )
    tag = head_tag()
    if tag != f"v{version}":
        fail(f"built version is {version} but HEAD is tagged {tag}")
    print(f"Release check passed: dist/ holds {version}, matching tag {tag}")


if __name__ == "__main__":
    main()
