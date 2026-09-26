"""End-to-end tests for the output writers and the step-skipping flags.

These runs are deliberately tiny: the goal is to exercise the PDB, HTML, and
SDF writers and the `skip_*` branches, not to check chemistry.
"""

import glob
import os

from rdkit import Chem

from gypsum_dl.start import prepare_molecules


def test_pdb_and_html_outputs_are_written(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\nc1ccccc1\tbenzene\n")
    output_folder = tmp_path / "out"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(output_folder),
            "job_manager": "serial",
            "separate_output_files": True,
            "add_pdb_output": True,
            "add_html_output": True,
            "max_variants_per_compound": 1,
            "thoroughness": 1,
            "use_durrant_lab_filters": True,
        }
    )
    pdb_files = glob.glob(os.path.join(str(output_folder), "*.pdb"))
    assert pdb_files
    with open(pdb_files[0]) as f:
        pdb_text = f.read()
    assert "REMARK Original SMILES string:" in pdb_text
    assert "REMARK Final SMILES string:" in pdb_text
    with open(os.path.join(str(output_folder), "gypsum_dl_success.html")) as f:
        assert "<div" in f.read()
    assert glob.glob(os.path.join(str(output_folder), "*.sdf"))


def test_skip_flags_produce_two_dimensional_output(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    output_folder = tmp_path / "out2d"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(output_folder),
            "job_manager": "serial",
            "num_processors": 4,
            "2d_output_only": True,
            "skip_adding_hydrogen": True,
            "skip_making_tautomers": True,
            "skip_enumerate_chiral_mol": True,
            "skip_enumerate_double_bonds": True,
            "skip_optimize_geometry": True,
            "skip_alternate_ring_conformations": True,
            "max_variants_per_compound": 1,
            "thoroughness": 1,
        }
    )
    assert os.path.exists(os.path.join(str(output_folder), "gypsum_dl_success.sdf"))


def test_unassigned_bond_and_unparseable_smiles_are_dropped(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("moosedogfacecat\tgarbage\nCCO\tethanol\n")
    output_folder = tmp_path / "out_mixed"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(output_folder),
            "job_manager": "serial",
            "2d_output_only": True,
            "max_variants_per_compound": 1,
            "thoroughness": 1,
        }
    )
    sdf_path = os.path.join(str(output_folder), "gypsum_dl_success.sdf")
    with open(sdf_path) as f:
        assert "garbage" not in f.read()


def test_nested_output_folder_is_created(tmp_path) -> None:
    # Regression (M5): os.mkdir raised FileNotFoundError when a parent segment
    # of output_folder was missing, before the "couldn't be created" message
    # could run. os.makedirs must build the whole tree.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    output_folder = tmp_path / "a" / "b" / "c"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(output_folder),
            "job_manager": "serial",
            "2d_output_only": True,
            "max_variants_per_compound": 1,
            "thoroughness": 1,
        }
    )
    assert os.path.isdir(str(output_folder))
    assert os.path.exists(os.path.join(str(output_folder), "gypsum_dl_success.sdf"))


def test_2d_output_has_nonzero_depiction_coordinates(tmp_path) -> None:
    # Regression (M6): with 2d_output_only, no conformer was ever loaded, so
    # SDWriter emitted a coordinate block of zeros. The SDF must carry real 2D
    # depiction coordinates, not all atoms stacked at the origin.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    output_folder = tmp_path / "out2dcoords"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(output_folder),
            "job_manager": "serial",
            "2d_output_only": True,
            "max_variants_per_compound": 1,
            "thoroughness": 1,
        }
    )
    sdf_path = os.path.join(str(output_folder), "gypsum_dl_success.sdf")
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    # The first SDF record is an empty placeholder holding the run parameters;
    # skip it and any other atomless record.
    mols = [m for m in supplier if m is not None and m.GetNumAtoms() > 0]
    assert mols
    conf = mols[0].GetConformer()
    assert any(
        abs(conf.GetAtomPosition(i).x) > 1e-6 or abs(conf.GetAtomPosition(i).y) > 1e-6
        for i in range(mols[0].GetNumAtoms())
    )