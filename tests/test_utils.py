"""Unit tests for the logging, sampling, and bookkeeping helpers in utils."""

import os
import subprocess
import sys

import pytest

from gypsum_dl import MyMol, utils
from gypsum_dl.MolContainer import MolContainer


def _mol(smiles: str, contnr_idx: int) -> MyMol.MyMol:
    """Build a MyMol tagged with a container index.

    The grouping helpers key off `contnr_idx`, which is normally set by the
    pipeline rather than the constructor.

    Args:
        smiles: SMILES string for the molecule.
        contnr_idx: Container index to attach.

    Returns:
        The tagged MyMol instance.
    """
    mol = MyMol.MyMol(smiles)
    mol.contnr_idx = contnr_idx
    return mol


def test_slug_empty() -> None:
    assert utils.slug("") == "untitled"


def test_slug_replaces_invalid_characters() -> None:
    assert utils.slug("ligand 1/2*") == "ligand_1_2_"


def test_slug_keeps_safe_characters() -> None:
    assert utils.slug("lig-and_1.2") == "lig-and_1.2"


def test_random_sample_truncates(capsys: pytest.CaptureFixture[str]) -> None:
    assert len(utils.random_sample([1, 2, 3, 4, 5], 2, "trimmed")) == 2
    assert "trimmed" in capsys.readouterr().out


def test_random_sample_keeps_everything_when_num_is_large() -> None:
    assert sorted(utils.random_sample([1, 2, 3], 10)) == [1, 2, 3]


def test_random_sample_tolerates_unhashable_items() -> None:
    assert len(utils.random_sample([[1], [2]], 5)) == 2


def test_group_mols_by_container_index_groups_and_skips_none() -> None:
    grouped = utils.group_mols_by_container_index(
        [_mol("CCO", 0), _mol("CCC", 1), None]
    )
    assert sorted(grouped.keys()) == [0, 1]
    assert len(grouped[0]) == 1


def test_group_mols_by_container_index_deduplicates() -> None:
    grouped = utils.group_mols_by_container_index([_mol("CCO", 0), _mol("OCC", 0)])
    assert len(grouped[0]) == 1


def test_fnd_contnrs_not_represntd_reports_gaps() -> None:
    contnrs = [MolContainer("CCO", "a", 0, {}), MolContainer("CCC", "b", 1, {})]
    assert utils.fnd_contnrs_not_represntd(contnrs, [_mol("CCO", 0)]) == [1]


def test_fnd_contnrs_not_represntd_empty_when_all_present() -> None:
    contnrs = [MolContainer("CCO", "a", 0, {})]
    assert utils.fnd_contnrs_not_represntd(contnrs, [_mol("CCO", 0)]) == []


def test_print_current_smiles_lists_container_contents(
    capsys: pytest.CaptureFixture[str],
) -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    utils.print_current_smiles([contnr])
    out = capsys.readouterr().out
    assert "MolContainer #0" in out
    assert "CCO" in out


def test_log_wraps_and_appends_trailing_whitespace(
    capsys: pytest.CaptureFixture[str],
) -> None:
    utils.log("\tindented message", trailing_whitespace="\n")
    out = capsys.readouterr().out
    assert "indented message" in out
    assert out.endswith("\n\n")


def test_log_preserves_embedded_newlines(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Regression: log() used a single textwrap.fill over the whole message,
    # which collapsed embedded newlines into one reflowed 80-column paragraph.
    # deal_with_failed_molecules logs "\n".join(failed_ones), so a long SMILES
    # list was mangled (SMILES broken mid-string, entries merged). Each input
    # line must survive as its own output line. (textwrap collapses runs of
    # whitespace, so tabs are not asserted on here.)
    smi_a = "C" * 50 + "O"
    smi_b = "C" * 50 + "N"
    utils.log("\n".join([smi_a, smi_b]))
    lines = capsys.readouterr().out.splitlines()
    assert smi_a in lines
    assert smi_b in lines


def test_exception_raises_with_message() -> None:
    with pytest.raises(Exception, match="boom"):
        utils.exception("boom")

def test_group_mols_by_container_index_dedups_in_first_seen_order() -> None:
    # Regression: deduplication used list(set(...)), whose order follows
    # PYTHONHASHSEED because MyMol hashes its canonical SMILES string. These
    # lists are the input to the seeded sampling in random_sample, so a seeded
    # run still selected different variants on every invocation.
    mols = [_mol(s, 0) for s in ["CCO", "CCC", "CCCC", "CCCCC", "OCC"]]

    grouped = utils.group_mols_by_container_index(mols)

    assert [m.smiles() for m in grouped[0]] == [m.smiles() for m in mols[:4]]


_HASH_SEED_SCRIPT = """
import random

from gypsum_dl import utils

random.seed(20260926)
print(",".join(utils.random_sample(["mol%d" % i for i in range(20)] * 2, 5)))
"""


def test_random_sample_ignores_the_interpreter_hash_seed() -> None:
    # Regression: random_sample deduplicated with set(), so the list handed to
    # random.shuffle was ordered differently in every interpreter invocation.
    # Which variants survived the max_variants_per_compound trim therefore
    # changed from run to run even with random_seed fixed. Run the same seeded
    # sample under several hash seeds; the selection must not move.
    selections = set()
    for hash_seed in ("1", "2", "3"):
        completed = subprocess.run(
            [sys.executable, "-c", _HASH_SEED_SCRIPT],
            capture_output=True,
            check=True,
            env=dict(os.environ, PYTHONHASHSEED=hash_seed),
            text=True,
        )
        selections.add(completed.stdout.strip().splitlines()[-1])

    assert len(selections) == 1