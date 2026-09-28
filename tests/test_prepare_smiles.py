"""Unit tests for the SMILES-preparation driver."""

from gypsum_dl import utils
from gypsum_dl.steps.smiles import PrepareSmiles


def _skip_all_params() -> dict:
    """Params that make prepare_smiles take every 'skip' branch.

    Keeps the driver from doing any real chemistry so only the debug-dump
    behavior is exercised.
    """
    return {
        "min_ph": 7.0,
        "max_ph": 7.0,
        "pka_precision": 1.0,
        "max_variants_per_compound": 1,
        "thoroughness": 1,
        "num_processors": 1,
        "job_manager": "serial",
        "let_tautomers_change_chirality": False,
        "Parallelizer": None,
        "skip_adding_hydrogen": True,
        "skip_making_tautomers": True,
        "use_durrant_lab_filters": False,
        "skip_enumerate_chiral_mol": True,
        "skip_enumerate_double_bonds": True,
    }


def test_prepare_smiles_debug_defaults_off(monkeypatch) -> None:
    # Regression: debug was hardcoded True, so every production job printed six
    # full print_current_smiles dumps. It now defaults to False via params.
    calls = []
    monkeypatch.setattr(PrepareSmiles, "desalt_orig_smi", lambda *a, **k: None)
    monkeypatch.setattr(utils, "print_current_smiles", lambda contnrs: calls.append(1))

    PrepareSmiles.prepare_smiles([], _skip_all_params())

    assert calls == []


def test_prepare_smiles_debug_dumps_when_enabled(monkeypatch) -> None:
    # With debug explicitly on, the six diagnostic dumps must still fire.
    calls = []
    monkeypatch.setattr(PrepareSmiles, "desalt_orig_smi", lambda *a, **k: None)
    monkeypatch.setattr(utils, "print_current_smiles", lambda contnrs: calls.append(1))

    params = _skip_all_params()
    params["debug"] = True
    PrepareSmiles.prepare_smiles([], params)

    assert len(calls) == 6
