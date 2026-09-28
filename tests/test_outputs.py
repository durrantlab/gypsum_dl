"""End-to-end tests for the output writers and the step-skipping flags.

These runs are deliberately tiny: the goal is to exercise the PDB, HTML, and
SDF writers and the `skip_*` branches, not to check chemistry.
"""

import glob
import os
import subprocess
import sys
import textwrap
from datetime import datetime

import pytest
from rdkit import Chem

import gypsum_dl
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.start import prepare_molecules
from gypsum_dl.steps.io import ProcessOutput
from gypsum_dl.steps.io.SaveToSDF import save_to_sdf
from gypsum_dl.steps.io.Web2DOutput import web_2d_output


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


def test_output_numbering_follows_input_position(tmp_path) -> None:
    # Regression (F5): the container index was a count of the records that
    # survived parsing rather than the record's position in the input file, so
    # a rejected line shifted every external number down by one. The filename
    # suffix, the UniqueID property, and the failed-molecule filename all come
    # from that index, so the molecule on line 2 was reported as input 1 and
    # cross-referencing output against input silently misattributed it.
    src = tmp_path / "input.smi"
    src.write_text("moosedogfacecat\tgarbage\nCCO\tethanol\n")
    output_folder = tmp_path / "out_numbering"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(output_folder),
            "job_manager": "serial",
            "separate_output_files": True,
            "2d_output_only": True,
            "max_variants_per_compound": 1,
            "thoroughness": 1,
        }
    )

    sdf_path = os.path.join(str(output_folder), "ethanol__input2.sdf")
    assert os.path.exists(sdf_path)
    assert not os.path.exists(os.path.join(str(output_folder), "ethanol__input1.sdf"))

    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    mols = [m for m in supplier if m is not None and m.GetNumAtoms() > 0]
    assert mols
    for mol in mols:
        assert mol.GetProp("UniqueID").startswith("2_")


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


def test_params_sdf_is_qualified_per_container(tmp_path) -> None:
    # Regression (M12): with separate_output_files=True (which MPI mode forces),
    # every task opened the same "gypsum_dl_params.sdf" for writing. Under MPI
    # all ranks run save_to_sdf concurrently against that one path, so the
    # provenance record interleaves or truncates. Each task holds exactly one
    # container, so the params file must be qualified by that container's
    # original index. Simulate two ranks writing into one folder and assert
    # they land in distinct files with nothing at the bare name.
    output_folder = tmp_path / "params_out"
    output_folder.mkdir()
    for idx in (0, 1):
        contnr = MolContainer("CCO", "ethanol", idx, {})
        contnr.add_smiles("CCO")
        save_to_sdf([contnr], {"thoroughness": 1}, True, str(output_folder))

    param_files = sorted(
        os.path.basename(p)
        for p in glob.glob(os.path.join(str(output_folder), "gypsum_dl_params*.sdf"))
    )
    assert param_files == [
        "gypsum_dl_params__input1.sdf",
        "gypsum_dl_params__input2.sdf",
    ]
    assert not os.path.exists(os.path.join(str(output_folder), "gypsum_dl_params.sdf"))


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


def test_genealogy_records_the_input_smiles(tmp_path) -> None:
    # Regression (bug 15): only the ionization-failure fallback wrote a
    # "(source)" entry, so the Genealogy field of a successfully processed
    # molecule began at "(protonated)" and never named the input SMILES the
    # variant derives from. The README points users at this field to trace a
    # problematic form back through the steps, so the source has to be there.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    output_folder = tmp_path / "out_genealogy"
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
    for mol in mols:
        assert mol.GetProp("Genealogy").startswith("CCO (source)")


def test_input_sdf_tags_do_not_replace_computed_output_fields(tmp_path) -> None:
    # Regression: add_container_properties merged the input record's tags over
    # mol_props at save time, after every step had already computed its own
    # values, and set_all_rdkit_mol_props wrote SMILES before iterating
    # mol_props. Feeding in an SDF that had been scored before (docking output,
    # a ChEMBL-style export) therefore reported the input's Energy and its
    # salted input SMILES as though Gypsum-DL had produced them.
    mol = Chem.MolFromSmiles("CCO.[Na+]")
    mol.SetProp("_Name", "lig")
    mol.SetDoubleProp("Energy", -9.5)
    mol.SetProp("SMILES", "CCO.[Na+]")
    mol.SetProp("UniqueID", "from_input")
    src = tmp_path / "input.sdf"
    writer = Chem.SDWriter(str(src))
    writer.write(mol)
    writer.close()

    output_folder = tmp_path / "out_input_props"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(output_folder),
            "job_manager": "serial",
            "max_variants_per_compound": 1,
            "thoroughness": 1,
        }
    )

    mols = _molecules_in(os.path.join(str(output_folder), "gypsum_dl_success.sdf"))
    assert mols
    for out in mols:
        # The variant is desalted, so its own SMILES carries no sodium.
        assert "Na" not in out.GetProp("SMILES")
        assert Chem.CanonSmiles(out.GetProp("SMILES")) == Chem.CanonSmiles("CCO")
        assert float(out.GetProp("Energy")) != pytest.approx(-9.5)
        assert out.GetProp("UniqueID") == "1_1"


def test_web_2d_output_declares_utf8_and_round_trips(tmp_path) -> None:
    # Regression (bug 17): the HTML was written with the platform default
    # encoding and carried no charset declaration, so a non-ASCII ligand name
    # rendered as mojibake even where the write itself succeeded.
    contnr = MolContainer("CCO", "caf\u00e9", 0, {})
    contnr.add_smiles("CCO")

    web_2d_output([contnr], str(tmp_path))

    html = (tmp_path / "gypsum_dl_success.html").read_text(encoding="utf-8")
    assert html.startswith('<meta charset="utf-8">')
    assert "caf\u00e9" in html


def test_web_2d_output_survives_a_non_utf8_locale(tmp_path) -> None:
    # Regression (bug 17): without an explicit encoding, writing a ligand named
    # "cafe\u0301" raised UnicodeEncodeError under a non-UTF-8 locale, at the very
    # end of a long run. The locale only affects the default encoding of a
    # freshly started interpreter, so this has to run in a child process.
    script = textwrap.dedent(
        """
        import sys

        from gypsum_dl.MolContainer import MolContainer
        from gypsum_dl.steps.io.Web2DOutput import web_2d_output

        contnr = MolContainer("CCO", "caf\\u00e9", 0, {})
        contnr.add_smiles("CCO")
        web_2d_output([contnr], sys.argv[1])
        """
    )

    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(gypsum_dl.__file__)))
    env = dict(os.environ)
    env["LC_ALL"] = "C"
    env["LANG"] = "C"
    # Keep CPython from restoring a UTF-8 default under the C locale (PEP 538
    # coercion, PEP 540 UTF-8 mode), which would let the child pass regardless.
    env["PYTHONCOERCECLOCALE"] = "0"
    env["PYTHONUTF8"] = "0"
    env.pop("PYTHONIOENCODING", None)
    existing_path = env.get("PYTHONPATH", "")
    env["PYTHONPATH"] = (
        repo_root + os.pathsep + existing_path if existing_path else repo_root
    )

    result = subprocess.run(
        [sys.executable, "-c", script, str(tmp_path)],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr

    written = (tmp_path / "gypsum_dl_success.html").read_bytes()
    assert written.startswith(b'<meta charset="utf-8">')
    # The name has to land as UTF-8 bytes, not as whatever the locale implied.
    assert "caf\u00e9".encode() in written


def _molecules_in(sdf_path: str) -> list:
    """Read back the real molecules from an SDF written by save_to_sdf.

    The first record is always an atomless placeholder holding the run
    parameters, so tests that count output molecules have to drop it.

    Args:
        sdf_path: Path to the SDF file to read.

    Returns:
        The records that carry atoms.
    """
    supplier = Chem.SDMolSupplier(sdf_path, removeHs=False)
    return [m for m in supplier if m is not None and m.GetNumAtoms() > 0]


def test_sdf_is_written_even_when_the_html_pass_fails(tmp_path, monkeypatch) -> None:
    # Regression: with --add_html_output, the unguarded 2D depiction pass ran
    # before the SDF was written. One variant that survives the pipeline but
    # cannot be depicted turned a completed run into no output at all, even
    # though the HTML is documented as a debugging aid.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")

    def boom(*args, **kwargs):
        raise RuntimeError("cannot depict this molecule")

    monkeypatch.setattr(ProcessOutput, "web_2d_output", boom)

    ProcessOutput.proccess_output(
        [contnr],
        {
            "separate_output_files": False,
            "output_folder": str(tmp_path),
            "add_pdb_output": False,
            "add_html_output": True,
        },
    )

    assert len(_molecules_in(str(tmp_path / "gypsum_dl_success.sdf"))) == 1


def test_save_to_sdf_skips_molecules_with_no_rdkit_mol(tmp_path) -> None:
    # Regression: SDWriter.write(None) raises, so one unwritable variant took
    # the whole SDF with it. Such a molecule has nothing to write anyway.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    contnr.add_smiles("CCCO")
    contnr.mols[0].rdkit_mol = None

    save_to_sdf([contnr], {"thoroughness": 1}, False, str(tmp_path))

    mols = _molecules_in(str(tmp_path / "gypsum_dl_success.sdf"))
    assert len(mols) == 1
    assert Chem.MolToSmiles(mols[0]) == "CCCO"


def test_save_to_sdf_finishes_the_file_when_a_molecule_raises(tmp_path) -> None:
    # Regression: the writer was flushed and closed only on the success path, so
    # an exception part way through the molecules left a truncated or empty SDF
    # behind. Whatever was written before the failure has to be readable.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    contnr.add_smiles("CCCO")

    def boom() -> None:
        raise RuntimeError("cannot load conformers")

    contnr.mols[1].load_conformers_into_rdkit_mol = boom

    with pytest.raises(RuntimeError):
        save_to_sdf([contnr], {"thoroughness": 1}, False, str(tmp_path))

    assert len(_molecules_in(str(tmp_path / "gypsum_dl_success.sdf"))) == 1


def test_zero_max_variants_still_writes_one_model_per_input(tmp_path) -> None:
    # Regression: max_variants_per_compound == 0 passes validation, and the
    # SMILES enumeration steps read it as "do not enumerate variants." The
    # ionization step had no such guard and the 3D steps read it as a cap of
    # zero survivors, so every container was emptied: the success SDF held
    # nothing but the parameters record, every input landed in
    # gypsum_dl_failed.smi, and the log blamed conformer generation.
    # Cyclohexane is here to carry the run through the non-aromatic
    # ring-conformer step as well.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\nC1CCCCC1\tcyclohexane\n")
    output_folder = tmp_path / "out_zero_variants"
    prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(output_folder),
            "job_manager": "serial",
            "max_variants_per_compound": 0,
            "thoroughness": 1,
        }
    )

    mols = _molecules_in(os.path.join(str(output_folder), "gypsum_dl_success.sdf"))
    assert len(mols) == 2
    assert all(m.GetNumAtoms() > 0 for m in mols)
    assert not os.path.exists(os.path.join(str(output_folder), "gypsum_dl_failed.smi"))


def _params_record_text(sdf_path: str) -> str:
    """Return the text of the parameters record an output SDF opens with.

    That record is an atomless placeholder, so reading it back through
    SDMolSupplier would depend on how RDKit treats a molecule with no atoms.
    Slicing the file text at the first record terminator avoids the question.

    Args:
        sdf_path: Path to the SDF file to read.

    Returns:
        The text of the first record, without its "$$$$" terminator.
    """
    with open(sdf_path, encoding="utf-8") as f:
        return f.read().split("$$$$")[0]


def test_params_record_carries_a_real_start_time(tmp_path) -> None:
    # Regression: the parameters record is written during the run, but
    # start_time, end_time, and run_time were assigned to params only after the
    # run returned, so every output SDF reported all three as their default 0.
    # The two that cannot be known when the record is written are no longer
    # written at all, and Parallelizer goes with them: it stringifies to an
    # address that changes every run, which defeats byte-identical reruns.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    output_folder = tmp_path / "out_params_record"
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

    record = _params_record_text(
        os.path.join(str(output_folder), "gypsum_dl_success.sdf")
    )
    assert "<run_time>" not in record
    assert "<end_time>" not in record
    assert "<Parallelizer>" not in record

    lines = record.splitlines()
    start_time_line = next(i for i, line in enumerate(lines) if "<start_time>" in line)
    # Raises if the recorded value is the "0" default rather than a timestamp.
    datetime.fromisoformat(lines[start_time_line + 1].strip())


def test_web_2d_output_skips_variants_the_other_writers_skip(tmp_path) -> None:
    # Regression: the HTML writer had none of the guards the SDF and PDB
    # writers have, so a variant with no rdkit_mol, or one whose no-hydrogen
    # SMILES could not be determined, raised part way through. The file is
    # opened with "w" and the caller swallows the exception, so the result was
    # a truncated HTML file that opens in a browser and looks complete.
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    contnr.add_smiles("CCCO")
    contnr.add_smiles("CCCCO")
    contnr.mols[1].rdkit_mol = None
    contnr.mols[2].can_smi_noh = None

    web_2d_output([contnr], str(tmp_path))

    html = (tmp_path / "gypsum_dl_success.html").read_text(encoding="utf-8")
    assert html.count('<div style="float: left') == 1
    assert "CCO" in html


def test_web_2d_output_escapes_the_ligand_name(tmp_path) -> None:
    # Regression (F8): the name was interpolated into the title attribute with
    # no escaping. Names come from an SDF _Name field or the tail of an SMI
    # line, so a quotation mark closed the attribute early and the rest of the
    # div was parsed as attributes, wrecking the depiction grid from that point
    # on.
    contnr = MolContainer("CCO", 'lig "A"', 0, {})
    contnr.add_smiles("CCO")

    web_2d_output([contnr], str(tmp_path))

    html_text = (tmp_path / "gypsum_dl_success.html").read_text(encoding="utf-8")
    assert 'title="lig &quot;A&quot;"' in html_text
    assert html_text.count('<div style="float: left') == 1
