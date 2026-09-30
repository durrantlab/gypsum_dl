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
from rdkit.ForceField.rdForceField import ForceField
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


def _embedded_pose() -> Chem.Mol:
    """Embed the target once, reproducibly, as the source of every reference.

    Returns:
        The embedded heavy-atom molecule.
    """
    full = Chem.AddHs(Chem.MolFromSmiles(TARGET))
    assert AllChem.EmbedMolecule(full, randomSeed=42) == 0
    return Chem.RemoveHs(full)


def _references(tmp_path: str, shift: float) -> list[str]:
    """Split one embedded pose of the target into two reference files.

    Args:
        tmp_path: Directory for the reference files.
        shift: Translation applied to the second reference only.

    Returns:
        Paths to the phenyl-side and pyridine-side references, in that order.
    """
    full = _embedded_pose()
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


def _write_probe_pdb(point: Point3D, path: str) -> str:
    """Write a one-atom receptor, so a clash can be placed exactly.

    Args:
        point: Where to put the probe carbon.
        path: Where to write the PDB file.

    Returns:
        The path written.
    """
    line = (
        "HETATM    1  C1  PRB A   1    "
        f"{point.x:8.3f}{point.y:8.3f}{point.z:8.3f}"
        "  1.00  0.00           C  \n"
    )
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(line + "END\n")
    return path


def test_receptor_pushes_a_clashing_hydrogen_away(tmp_path: Path) -> None:
    # Regression: CombineMols leaves ligand and receptor as separate
    # fragments, and UFF drops interfragment terms by default, so the
    # receptor had no effect on the minimized ligand.
    script = _load_script()
    refs = _references(str(tmp_path), 0.0)

    # Reference atom 0 is an aromatic CH (ring neighbors 1 and 5), so its
    # hydrogen lies on the exterior bisector about 1.08 A out. The probe sits
    # 0.8 A beyond that. Whichever target CH the MCS maps onto this position
    # is held fixed there, so the geometry does not depend on the mapping.
    ref_conf = Chem.MolFromMolFile(refs[0]).GetConformer()
    carbon = ref_conf.GetAtomPosition(0)
    outward = carbon * 2 - ref_conf.GetAtomPosition(1) - ref_conf.GetAtomPosition(5)
    outward.Normalize()
    probe = carbon + outward * (1.08 + 0.8)
    receptor = _write_probe_pdb(probe, os.path.join(tmp_path, "probe.pdb"))
    out = os.path.join(tmp_path, "out.sdf")

    script.generate_conformer(refs, TARGET, out, receptor)

    result = Chem.MolFromMolFile(out, removeHs=False)
    assert result is not None
    conf = result.GetConformer()
    fixed_carbon = [
        atom
        for atom in result.GetAtoms()
        if conf.GetAtomPosition(atom.GetIdx()).Distance(carbon) < 0.05
    ]
    assert len(fixed_carbon) == 1
    hydrogens = [n for n in fixed_carbon[0].GetNeighbors() if n.GetAtomicNum() == 1]
    assert len(hydrogens) == 1

    # Without receptor terms the hydrogen settles about 0.8 A from the probe.
    # Its carbon is fixed, so bending is all it can do, and against UFF's
    # angle terms that takes it to roughly 1.7 A rather than past 2 A.
    distance = conf.GetAtomPosition(hydrogens[0].GetIdx()).Distance(probe)
    assert distance > 1.3


class _FrozenForceField:
    """Wrap a force field so minimization leaves the coordinates untouched.

    That exposes the starting geometry the script hands to the minimizer,
    which is what has to be in the reference frame.
    """

    def __init__(self, ff: ForceField) -> None:
        self._ff = ff

    def AddFixedPoint(self, idx: int) -> None:
        """Pass fixed points through, as the script expects."""
        self._ff.AddFixedPoint(idx)

    def Minimize(self, maxIts: int = 200) -> int:
        """Skip minimization and report convergence."""
        return 0


def test_unmapped_atoms_start_in_the_reference_frame(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Regression: ETKDG embeds near the origin, and only the mapped atoms were
    # overwritten with reference coordinates. With a reference far from the
    # origin, every unmapped atom (all hydrogens, and here the whole pyridine)
    # started tens of A from the atoms it is bonded to.
    script = _load_script()
    ref = _write_fragment(
        _embedded_pose(), PHENYL_SIDE, 30.0, os.path.join(tmp_path, "a.mol")
    )
    original = script.AllChem.MMFFGetMoleculeForceField

    def frozen(*args: object, **kwargs: object) -> _FrozenForceField:
        """Build the real force field, then freeze it."""
        return _FrozenForceField(original(*args, **kwargs))

    monkeypatch.setattr(script.AllChem, "MMFFGetMoleculeForceField", frozen)
    out = os.path.join(tmp_path, "out.sdf")

    script.generate_conformer([ref], TARGET, out)

    result = Chem.MolFromMolFile(out, removeHs=False)
    assert result is not None
    conf = result.GetConformer()
    lengths = [
        conf.GetAtomPosition(b.GetBeginAtomIdx()).Distance(
            conf.GetAtomPosition(b.GetEndAtomIdx())
        )
        for b in result.GetBonds()
    ]
    # Unminimized, so bonds next to mapped atoms carry the alignment residual
    # between two conformers (the linker torsion differs). That stays well
    # under 4 A, while the frame offset put these bonds near 30 A.
    assert max(lengths) < 4.0
