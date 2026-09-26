"""End-to-end tests for the output writers and the step-skipping flags.

These runs are deliberately tiny: the goal is to exercise the PDB, HTML, and
SDF writers and the `skip_*` branches, not to check chemistry.
"""

import glob
import os

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