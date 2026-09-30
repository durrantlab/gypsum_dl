"""Regression tests for restrained_min_tests/test.py.

The script copies each reference's coordinates as-is. References in different
frames therefore put their fragments far apart, and because every copied atom
is held fixed during minimization, the joining bond stayed stretched in an
output written without any error.
"""

import importlib.util
import os
from pathlib import Path
from types import ModuleType

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

SCRIPT = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "restrained_min_tests", "test.py"
)

# Atoms 0-5 are the phenyl ring, 6-7 the ethylene linker, 8-13 the pyridine.
TARGET = "c1ccc(cc1)CCc1ccncc1"
PHENYL_SIDE = list(range(8))
PYRIDINE_SIDE = list(range(8, 14))
JUNCTION = (7, 8)


def _load_script() -> ModuleType:
    """Import the script by path, since its directory is not a package.

    Returns:
        The loaded module.
    """
    spec = importlib.util.spec_from_file_location("restrained_min_script", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_fragment(full: Chem.Mol, keep: list[int], shift: float, path: str) -> str:
    """Write part of an embedded molecule as a reference, optionally moved.

    Args:
        full: Embedded heavy-atom molecule to cut the fragment from.
        keep: Atom indices that make up the fragment.
        shift: Translation along x, in angstroms, applied to the fragment.
        path: Where to write the MOL file.

    Returns:
        The path written.
    """
    frag = Chem.RWMol(full)
    for idx in sorted(set(range(full.GetNumAtoms())) - set(keep), reverse=True):
        frag.RemoveAtom(idx)
    frag = frag.GetMol()
    Chem.SanitizeMol(frag)
    conf = frag.GetConformer()
    for idx in range(frag.GetNumAtoms()):
        pos = conf.GetAtomPosition(idx)
        conf.SetAtomPosition(idx, Point3D(pos.x + shift, pos.y, pos.z))
    Chem.MolToMolFile(frag, path)
    return path


def _references(tmp_path: str, shift: float) -> list[str]:
    """Split one embedded pose of the target into two reference files.

    Args:
        tmp_path: Directory for the reference files.
        shift: Translation applied to the second reference only.

    Returns:
        Paths to the phenyl-side and pyridine-side references, in that order.
    """
    full = Chem.AddHs(Chem.MolFromSmiles(TARGET))
    assert AllChem.EmbedMolecule(full, randomSeed=42) == 0
    full = Chem.RemoveHs(full)
    return [
        _write_fragment(full, PHENYL_SIDE, 0.0, os.path.join(tmp_path, "a.mol")),
        _write_fragment(full, PYRIDINE_SIDE, shift, os.path.join(tmp_path, "b.mol")),
    ]


def test_references_in_one_frame_are_joined(tmp_path: Path) -> None:
    script = _load_script()
    out = os.path.join(tmp_path, "out.sdf")

    script.generate_conformer(_references(str(tmp_path), 0.0), TARGET, out)

    result = Chem.MolFromMolFile(out, removeHs=False)
    assert result is not None
    conf = result.GetConformer()
    length = conf.GetAtomPosition(JUNCTION[0]).Distance(
        conf.GetAtomPosition(JUNCTION[1])
    )
    assert length < 2.0


def test_references_in_different_frames_are_rejected(tmp_path: Path) -> None:
    script = _load_script()
    out = os.path.join(tmp_path, "out.sdf")

    with pytest.raises(SystemExit):
        script.generate_conformer(_references(str(tmp_path), 20.0), TARGET, out)

    assert not os.path.exists(out)


def test_find_cross_reference_gaps_ignores_bonds_within_one_reference() -> None:
    script = _load_script()
    mol = Chem.MolFromSmiles("CCC")
    positions = {0: Point3D(0, 0, 0), 1: Point3D(20, 0, 0), 2: Point3D(21.5, 0, 0)}

    # A long bond inside one reference is that reference's own geometry.
    assert script.find_cross_reference_gaps(mol, positions, {0: 0, 1: 0, 2: 1}) == []

    gaps = script.find_cross_reference_gaps(mol, positions, {0: 0, 1: 1, 2: 1})
    assert [(a, b) for a, b, _ in gaps] == [(0, 1)]
