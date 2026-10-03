"""End-to-end checks that output stereochemistry is self-consistent and kept.

Unit tests exercise each enumeration step on molecules built one way, but the
pipeline builds them other ways (the double-bond fix first passed its unit
test while real runs still flipped a specified bond). These tests run the
whole pipeline once and inspect the SDF it writes.
"""

from collections import defaultdict

import pytest
from rdkit import Chem

from gypsum_dl.start import prepare_molecules

# Name to input SMILES. Covers specified and unspecified cis/trans bonds, a
# partially specified conjugated diene, specified and unspecified chiral
# centers (one that ionizes), ring stereo, and both kinds together.
STEREO_INPUTS: dict[str, str] = {
    "trans_butene": "C/C=C/C",
    "cis_butene": "C/C=C\\C",
    "trans_fluorochloroethene": "F/C=C/Cl",
    "unspecified_butene": "CC=CC",
    "partially_specified_diene": "C/C=C/C=CC",
    "l_alanine": "C[C@H](N)C(=O)O",
    "r_butanol": "C[C@@H](O)CC",
    "unspecified_butanol": "CC(O)CC",
    "methylcyclohexanol": "C[C@H]1CC[C@H](O)CC1",
    "alkene_and_center": "C/C=C/[C@H](C)O",
}


@pytest.fixture(scope="module")
def output_records(tmp_path_factory: pytest.TempPathFactory) -> list[Chem.Mol]:
    """Run the pipeline once on the stereo inputs and read every record.

    Module-scoped because the run is the expensive part, and every test here
    only reads its output.

    Args:
        tmp_path_factory: Pytest's factory for a module-lifetime directory.

    Returns:
        Each output record, read with its hydrogens so stereo can be
            perceived from the coordinates.
    """
    tmp = tmp_path_factory.mktemp("stereo")
    src = tmp / "stereo.smi"
    src.write_text("".join(f"{smi}\t{name}\n" for name, smi in STEREO_INPUTS.items()))
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(tmp),
            "separate_output_files": False,
            "job_manager": "serial",
            "max_variants_per_compound": 8,
            "thoroughness": 3,
            "use_durrant_lab_filters": True,
        }
    )
    return [
        m
        for m in Chem.SDMolSupplier(str(tmp / "gypsum_dl_success.sdf"), removeHs=False)
        if m is not None and m.HasProp("SMILES")
    ]


def _tags_by_name(records: list[Chem.Mol]) -> dict[str, set[str]]:
    """Group the canonical SMILES tags of the records by input name.

    Args:
        records: Output records.

    Returns:
        Input name to the set of canonical SMILES written for it.
    """
    tags: dict[str, set[str]] = defaultdict(set)
    for m in records:
        tags[m.GetProp("_Name")].add(Chem.CanonSmiles(m.GetProp("SMILES")))
    return tags


def test_every_input_produces_output(output_records: list[Chem.Mol]) -> None:
    """Guard the other tests against passing because records went missing."""
    assert set(_tags_by_name(output_records)) == set(STEREO_INPUTS)


def test_output_geometry_matches_the_smiles_tag(
    output_records: list[Chem.Mol],
) -> None:
    """Stereo perceived from each record's 3D coordinates must match its tag.

    The tag and the coordinates come from differently built RDKit molecules,
    so nothing guaranteed they agreed, and nothing tested it.
    """
    mismatches = []
    for m in output_records:
        perceived = Chem.Mol(m)
        Chem.AssignStereochemistryFrom3D(perceived)
        from_3d = Chem.MolToSmiles(Chem.RemoveHs(perceived))
        tag = Chem.CanonSmiles(m.GetProp("SMILES"))
        if from_3d != tag:
            mismatches.append((m.GetProp("_Name"), tag, from_3d))

    assert not mismatches


@pytest.mark.parametrize(
    "name", ["trans_butene", "cis_butene", "trans_fluorochloroethene", "r_butanol"]
)
def test_fully_specified_inputs_keep_their_stereo(
    output_records: list[Chem.Mol], name: str
) -> None:
    """An input with all of its stereo specified must come out unchanged.

    These inputs have no ionizable groups or tautomers, so the only output
    is the input itself.
    """
    assert _tags_by_name(output_records)[name] == {
        Chem.CanonSmiles(STEREO_INPUTS[name])
    }


def test_partially_specified_diene_keeps_its_specified_bond(
    output_records: list[Chem.Mol],
) -> None:
    """Only the unspecified bond of C/C=C/C=CC may be varied.

    Regression: the pipeline also wrote the (Z,Z) isomer, flipping the bond
    the input fixed as trans, because the single bond the two double bonds
    share was given a direction for the unspecified one.
    """
    assert _tags_by_name(output_records)["partially_specified_diene"] == {
        Chem.CanonSmiles("C/C=C/C=C/C"),
        Chem.CanonSmiles("C/C=C/C=C\\C"),
    }


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("unspecified_butene", {"C/C=C/C", "C/C=C\\C"}),
        ("unspecified_butanol", {"CC[C@H](C)O", "CC[C@@H](C)O"}),
    ],
)
def test_unspecified_stereo_is_enumerated(
    output_records: list[Chem.Mol], name: str, expected: set[str]
) -> None:
    """Each unspecified element must come out in both configurations."""
    assert _tags_by_name(output_records)[name] == {
        Chem.CanonSmiles(s) for s in expected
    }
