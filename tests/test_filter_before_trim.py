"""Regression tests for applying the Durrant-lab filters before variant trims.

The filter step used to see only the variants that survived the tautomer
step's trim to max_variants_per_compound. When that trim kept nothing but
forms the filters reject, the compound failed although an acceptable tautomer
had been generated and discarded. NSC 136306, an acylurea, is the reported
case: at pH 7.4 its one ionization state has a keto form the filters accept
and several enol and iminol forms they reject.
"""

import os
from pathlib import Path

import pytest
from rdkit import Chem

from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.parallelizer import seed_generators
from gypsum_dl.start import prepare_molecules
from gypsum_dl.steps.smiles import MakeTautomers
from gypsum_dl.steps.smiles.DurrantLabFilter import passes_durrant_lab_filters
from gypsum_dl.steps.smiles.MakeTautomers import make_tauts

ACYLUREA_INPUT = "O=C(NC(=O)C(Cl)Cl)Nc1cccc2c1-c1ccccc1C2"
ACYLUREA_KETO_ANION = "O=C([N-]C(=O)C(Cl)Cl)Nc1cccc2c1-c1ccccc1C2"
ACYLUREA_NAME = "136306"


def _canon(smiles: str) -> str:
    """Canonicalize a SMILES so written and expected forms compare equal.

    Args:
        smiles: Any valid SMILES.

    Returns:
        RDKit's canonical SMILES for it.
    """
    return Chem.CanonSmiles(smiles)


def _acylurea_anion_contnr(idx: int = 0) -> MolContainer:
    """Build the reported compound as the tautomer step would receive it.

    Ionization runs first, so the tautomer step sees the keto anion rather
    than the neutral input. Starting from it keeps these tests independent of
    Dimorphite-DL.

    Args:
        idx: Container index to assign.

    Returns:
        A container whose one variant is the keto anion.
    """
    contnr = MolContainer(ACYLUREA_INPUT, ACYLUREA_NAME, idx, {})
    contnr.add_smiles(ACYLUREA_KETO_ANION)
    return contnr


def _run_2d(tmp_path: Path, smi_line: str, use_filters: bool) -> Path:
    """Run the reported 2D configuration on one compound.

    Args:
        tmp_path: Pytest's per-test directory.
        smi_line: One line of a .smi file (SMILES, whitespace, name).
        use_filters: Whether to turn on the Durrant-lab filters.

    Returns:
        The output folder.
    """
    src = tmp_path / "input.smi"
    src.write_text(smi_line + "\n")
    output_folder = tmp_path / "out"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(output_folder),
            "job_manager": "serial",
            "2d_output_only": True,
            "max_variants_per_compound": 1,
            "min_ph": 7.4,
            "max_ph": 7.4,
            "pka_precision": 0,
            "use_durrant_lab_filters": use_filters,
        }
    )
    return output_folder


def _written_smiles(sdf_path: Path) -> list[str]:
    """Read the canonical SMILES of every compound record in an output SDF.

    Args:
        sdf_path: A Gypsum-DL SDF.

    Returns:
        One canonical SMILES per compound record, in file order. The
            parameter record, which has no SMILES property, is skipped.
    """
    return [
        _canon(m.GetProp("SMILES"))
        for m in Chem.SDMolSupplier(str(sdf_path))
        if m is not None and m.HasProp("SMILES")
    ]


def _failed_text(output_folder: Path) -> str:
    """Collect everything written to the failed-molecule files.

    Args:
        output_folder: A Gypsum-DL output folder.

    Returns:
        The concatenated contents of every gypsum_dl_failed*.smi file, or an
            empty string when none was written.
    """
    return "".join(
        p.read_text() for p in sorted(output_folder.glob("gypsum_dl_failed*.smi"))
    )


def test_acylurea_keeps_its_acceptable_tautomer_in_2d(tmp_path: Path) -> None:
    out = _run_2d(tmp_path, f"{ACYLUREA_INPUT} {ACYLUREA_NAME}", use_filters=True)

    assert _written_smiles(out / "gypsum_dl_success.sdf") == [
        _canon(ACYLUREA_KETO_ANION)
    ]
    assert ACYLUREA_NAME not in _failed_text(out)


def test_acylurea_produces_an_sdf_in_3d(tmp_path: Path) -> None:
    src = tmp_path / "input.smi"
    src.write_text(f"{ACYLUREA_INPUT} {ACYLUREA_NAME}\n")
    out = tmp_path / "out3d"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(out),
            "job_manager": "serial",
            "separate_output_files": True,
            "max_variants_per_compound": 1,
            "min_ph": 7.4,
            "max_ph": 7.4,
            "pka_precision": 0,
            "use_durrant_lab_filters": True,
        }
    )

    assert os.path.exists(out / f"{ACYLUREA_NAME}__input1.sdf")
    assert _failed_text(out) == ""


def test_compound_with_no_acceptable_variant_still_fails(tmp_path: Path) -> None:
    # Guards against a fix that quietly disables the filters: boron matches a
    # prohibited pattern in every form, so the compound must still fail, even
    # though the tautomer step's carry-over hands its structure forward.
    out = _run_2d(tmp_path, "B(O)(O)O boric_acid", use_filters=True)

    assert "boric_acid" in _failed_text(out)
    assert _written_smiles(out / "gypsum_dl_success.sdf") == []


@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_make_tauts_trims_to_an_acceptable_tautomer(seed: int) -> None:
    # The trim samples thoroughness * max_variants_per_compound candidates and
    # ranks them by energy, so several seeds are tried: before the filters ran
    # ahead of the trim, whether the keto form survived depended on both.
    seed_generators(seed)
    contnr = _acylurea_anion_contnr()

    make_tauts([contnr], 1, 3, 1, "serial", False, None, True)

    assert [_canon(m.smiles()) for m in contnr.mols] == [_canon(ACYLUREA_KETO_ANION)]


def test_make_tauts_without_filters_never_consults_them(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # With the filters off, the tautomer step must behave as it did before
    # they could run inside it: the predicate is never called, and rejected
    # forms remain eligible.
    def forbidden(*args: object) -> bool:
        """Fail the test if the tautomer step consults the filters.

        Args:
            *args: Whatever the step would have passed the real predicate.

        Returns:
            Never returns.
        """
        raise AssertionError("filters consulted with use_durrant_lab_filters off")

    monkeypatch.setattr(MakeTautomers, "passes_durrant_lab_filters", forbidden)
    contnr = _acylurea_anion_contnr()

    make_tauts([contnr], 8, 1, 1, "serial", False, None)

    assert any(not passes_durrant_lab_filters(m, ACYLUREA_INPUT) for m in contnr.mols)


# Ordinary compounds none of whose tautomers the filters reject.
UNAFFECTED_SPECS = [
    ("CC(=O)CC", "butanone"),
    ("CC(=O)CC(C)=O", "acetylacetone"),
    ("Oc1ccccc1", "phenol"),
    ("CCOC(C)=O", "ethyl_acetate"),
]


def _unaffected_contnrs() -> list[MolContainer]:
    """Build fresh, distinctly indexed containers for UNAFFECTED_SPECS.

    Each run needs its own containers, because make_tauts rewrites them and
    caches ranking energies on them.

    Returns:
        One single-variant container per spec, in spec order.
    """
    contnrs = []
    for idx, (smiles, name) in enumerate(UNAFFECTED_SPECS):
        contnr = MolContainer(smiles, name, idx, {})
        contnr.add_smiles(smiles)
        contnrs.append(contnr)
    return contnrs


def test_filters_do_not_change_compounds_with_no_rejected_tautomer() -> None:
    # For a compound none of whose tautomers the filters reject, the candidate
    # list is the same with the filters on or off, so a seeded run must keep
    # the same variants either way. Several containers go through one call so
    # that the shared random stream is compared too.
    wide = _unaffected_contnrs()
    make_tauts(wide, 50, 1, 1, "serial", False, None)
    for contnr in wide:
        assert all(passes_durrant_lab_filters(m, contnr.orig_smi) for m in contnr.mols)

    kept: dict[bool, list[list[str]]] = {}
    for filters in (False, True):
        seed_generators(7)
        contnrs = _unaffected_contnrs()
        make_tauts(contnrs, 1, 1, 1, "serial", False, None, filters)
        kept[filters] = [[m.smiles() for m in c.mols] for c in contnrs]

    assert kept[True] == kept[False]
