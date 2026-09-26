"""Unit tests for parameter validation and the small helpers in start.py."""

import json
import os
import random
from collections.abc import Callable

import numpy
import pytest

from gypsum_dl import start
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.steps.smiles.DeSaltOrigSmiles import desalt_orig_smi
from gypsum_dl.steps.smiles.DurrantLabFilter import durrant_lab_filters


def test_detect_unassigned_bonds_accepts_valid_smiles() -> None:
    assert start.detect_unassigned_bonds("CCO") == "CCO"


def test_detect_unassigned_bonds_rejects_garbage() -> None:
    assert start.detect_unassigned_bonds("moosedogfacecat") is None


def test_make_type_dict_maps_scalar_types() -> None:
    type_dict = start.make_type_dict({"a": 1, "b": 1.5, "c": True, "d": "x"})
    assert type_dict == {"a": int, "b": float, "c": bool, "d": str}


def test_make_type_dict_rejects_unsupported_types() -> None:
    with pytest.raises(Exception, match="No value can have type"):
        start.make_type_dict({"a": None})


def test_merge_parameters_promotes_int_to_float() -> None:
    default = {"min_ph": 6.4}
    start.merge_parameters(default, {"min_ph": 7})
    assert isinstance(default["min_ph"], float)
    assert default["min_ph"] == pytest.approx(7.0)


def test_merge_parameters_rejects_unknown_parameter() -> None:
    with pytest.raises(Exception, match="Unrecognized parameter"):
        start.merge_parameters({"min_ph": 6.4}, {"bogus": 1})


def test_merge_parameters_rejects_wrong_type() -> None:
    with pytest.raises(Exception, match="must be of"):
        start.merge_parameters({"min_ph": 6.4}, {"min_ph": "high"})


def test_set_parameters_lowercases_keys_and_fills_defaults(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params = start.set_parameters({"SOURCE": str(src), "job_manager": "SERIAL"})
    assert params["source"] == os.path.abspath(str(src))
    assert params["job_manager"] == "serial"
    assert params["thoroughness"] == 3
    assert params["max_variants_per_compound"] == 5
    # Regression (M10): a random_seed parameter must exist so runs can be made
    # reproducible; the default (-1) leaves the RNG unseeded.
    assert params["random_seed"] == -1


def test_set_parameters_accepts_random_seed(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params = start.set_parameters({"source": str(src), "random_seed": 42})
    assert params["random_seed"] == 42


def test_seed_random_number_generators_seeds_both_generators() -> None:
    # Regression: only the random module was seeded, so the numpy stream that
    # scipy's kmeans2 draws from stayed unseeded and a run given a random_seed
    # still varied between invocations.
    start.seed_random_number_generators({"random_seed": 42})
    first = (random.random(), float(numpy.random.random()))

    start.seed_random_number_generators({"random_seed": 42})
    second = (random.random(), float(numpy.random.random()))

    assert first == second


def test_seed_random_number_generators_leaves_generators_alone_when_negative() -> None:
    random.seed(7)
    numpy.random.seed(7)
    expected = (random.random(), float(numpy.random.random()))

    random.seed(7)
    numpy.random.seed(7)
    start.seed_random_number_generators({"random_seed": -1})

    assert (random.random(), float(numpy.random.random())) == expected


def test_seed_random_number_generators_accepts_seed_numpy_cannot_take() -> None:
    # numpy only accepts seeds that fit in 32 bits, so a larger one must not
    # take down the run.
    start.seed_random_number_generators({"random_seed": 2**40 + 1})


def test_set_parameters_rejects_zero_thoroughness(tmp_path) -> None:
    # Regression: only the type of thoroughness was checked. A zero reached
    # math.log(thoroughness * max_variants_per_compound, 2) inside a worker,
    # where the ValueError was swallowed and the molecule silently dropped.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    with pytest.raises(Exception, match="thoroughness"):
        start.set_parameters({"source": str(src), "thoroughness": 0})


def test_set_parameters_rejects_negative_thoroughness(tmp_path) -> None:
    # Regression: a negative thoroughness ran to completion. random_sample
    # slices with lst[:num] and pick_lowest_enrgy_mols multiplies by
    # thoroughness, so the candidate pools were quietly truncated and the run
    # wrote a full-looking SDF with no warning anywhere.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    with pytest.raises(Exception, match="thoroughness"):
        start.set_parameters({"source": str(src), "thoroughness": -1})


def test_set_parameters_rejects_negative_max_variants(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    with pytest.raises(Exception, match="max_variants_per_compound"):
        start.set_parameters({"source": str(src), "max_variants_per_compound": -1})


def test_set_parameters_allows_zero_max_variants(tmp_path) -> None:
    # Zero is the sentinel the SMILES enumeration steps check to skip
    # themselves, so validation must not reject it.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params = start.set_parameters({"source": str(src), "max_variants_per_compound": 0})
    assert params["max_variants_per_compound"] == 0


def test_set_parameters_allows_thoroughness_of_one(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params = start.set_parameters({"source": str(src), "thoroughness": 1})
    assert params["thoroughness"] == 1


def test_finalize_params_requires_source() -> None:
    with pytest.raises(Exception, match="source"):
        start.finalize_params({"source": "", "output_folder": "./"})


def test_finalize_params_defaults_output_folder_next_to_source(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params = start.finalize_params(
        {
            "source": str(src),
            "output_folder": "",
            "add_pdb_output": False,
            "separate_output_files": False,
            "job_manager": "Serial",
        }
    )
    assert params["output_folder"].endswith(f"output{os.sep}")
    assert params["job_manager"] == "serial"


def test_finalize_params_derives_source_dir_by_dirname(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression (M11): source_dir used str.strip(basename), which strips a
    # character *set* off both ends rather than removing a suffix. On POSIX
    # absolute paths the leading "/" and the last "/" hide the bug, but any
    # path whose leading directory character also appears in the basename gets
    # mangled (e.g. "smiles_dir/mol.smi" -> "es_dir/", and on Windows the drive
    # letter is eaten). Pin abspath to identity so a relative path reaches the
    # source_dir computation intact, then assert it equals os.path.dirname.
    # isfile is pinned too, because finalize_params now rejects a source that
    # does not exist and this synthetic path never does.
    monkeypatch.setattr(start.os.path, "abspath", lambda p: p)
    monkeypatch.setattr(start.os.path, "isfile", lambda p: True)
    params = start.finalize_params(
        {
            "source": "smiles_dir/mol.smi",
            "output_folder": "",
            "add_pdb_output": False,
            "separate_output_files": False,
            "job_manager": "serial",
        }
    )
    assert params["output_folder"] == "smiles_dir" + os.sep + "output" + os.sep


def test_add_mol_id_props_assigns_unique_ids() -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    start.add_mol_id_props([contnr])
    assert contnr.mols[0].rdkit_mol.GetProp("UniqueID") == "1_1"


def test_add_mol_id_props_ids_differ_across_separate_calls() -> None:
    # Regression (bug 16): the id came from a counter that restarted at 1 in
    # every call, and execute_gypsum_dl (which calls this) runs once per mpi
    # rank and once per task in separate_output_files mode. Concatenating the
    # per-input SDFs therefore produced many records labeled "1". Simulate two
    # such tasks and assert their ids do not collide.
    first = MolContainer("CCO", "ethanol", 0, {})
    first.add_smiles("CCO")
    second = MolContainer("CCCO", "propanol", 1, {})
    second.add_smiles("CCCO")

    start.add_mol_id_props([first])
    start.add_mol_id_props([second])

    first_id = first.mols[0].rdkit_mol.GetProp("UniqueID")
    second_id = second.mols[0].rdkit_mol.GetProp("UniqueID")
    assert first_id != second_id


def test_add_mol_id_props_ids_survive_mpi_renumbering() -> None:
    # Regression (bug 16): mpi mode renumbers each container to index zero of
    # its own job (contnr.update_idx(0)), so the id has to be built from
    # contnr_idx_orig rather than the working index.
    contnr = MolContainer("CCO", "ethanol", 3, {})
    contnr.add_smiles("CCO")
    contnr.update_idx(0)

    start.add_mol_id_props([contnr])

    assert contnr.mols[0].rdkit_mol.GetProp("UniqueID") == "4_1"


def test_deal_with_failed_molecules_writes_failure_file(tmp_path) -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    start.deal_with_failed_molecules([contnr], {"output_folder": str(tmp_path)})
    assert "ethanol" in (tmp_path / "gypsum_dl_failed.smi").read_text()


def test_deal_with_failed_molecules_skips_when_nothing_failed(tmp_path) -> None:
    contnr = MolContainer("CCO", "ethanol", 0, {})
    contnr.add_smiles("CCO")
    start.deal_with_failed_molecules([contnr], {"output_folder": str(tmp_path)})
    assert not os.path.exists(os.path.join(str(tmp_path), "gypsum_dl_failed.smi"))


def test_prepare_molecules_rejects_malformed_json(tmp_path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text("{not json")
    with pytest.raises(Exception, match="properly formed"):
        start.prepare_molecules({"json": str(bad)})


def test_prepare_molecules_rejects_missing_json(tmp_path) -> None:
    # Regression (B17): the JSON was read with json.load(open(...)) inside a
    # bare `except:`. The bare except is now `except (OSError, ValueError)`, so
    # a missing file still surfaces the friendly error rather than being caught
    # by an over-broad handler (and the file handle no longer leaks).
    with pytest.raises(Exception, match="properly formed"):
        start.prepare_molecules({"json": str(tmp_path / "does_not_exist.json")})


def test_set_parameters_num_processors_defaults_to_all_cores(tmp_path) -> None:
    # Regression (B18): set_parameters defaults num_processors to -1 (all
    # cores). The CLI must not override this with a hardcoded 1.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params = start.set_parameters({"source": str(src)})
    assert params["num_processors"] == -1


def test_deal_with_failed_molecules_separates_files_per_container(tmp_path) -> None:
    # Regression: the failure list always went to one fixed filename, opened
    # for writing. In separate-file mode (which mpi mode forces) every task
    # runs this, so each one truncated the previous task's report and only the
    # last writer's failures survived.
    for idx, smiles, name in ((0, "CCO", "ethanol"), (1, "CCCO", "propanol")):
        contnr = MolContainer(smiles, name, idx, {})
        start.deal_with_failed_molecules(
            [contnr],
            {"output_folder": str(tmp_path), "separate_output_files": True},
        )

    assert "ethanol" in (tmp_path / "gypsum_dl_failed__input1.smi").read_text()
    assert "propanol" in (tmp_path / "gypsum_dl_failed__input2.smi").read_text()


def test_prepare_molecules_mpi_reindex_restamps_original_mol(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Regression: the mpi branch renumbers every container to zero, because
    # each one is then run in isolation. It did so by assigning contnr_idx
    # directly, which left mol_orig_frm_inp_smi.contnr_idx at the pre-mpi
    # value. The desalter returns that very object for single-fragment inputs,
    # so the stale index rode into contnr.mols and the Durrant-lab filter, which
    # regroups molecules by mol.contnr_idx, found nothing under the container's
    # key and emptied it.
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\nCCCO\tpropanol\n")

    captured: list[tuple[list[MolContainer], dict[str, object]]] = []

    class _StubParallelizer:
        """Stand in for Parallelizer so the mpi grouping branch can be reached.

        Real mpi mode needs mpi4py and an mpi launcher, but the branch that
        renumbers the containers is selected purely on return_mode(), so a stub
        that claims mpi and records its jobs exercises it in process.
        """

        def __init__(self, *args: object, **kwargs: object) -> None:
            pass

        def return_mode(self) -> str:
            """Report mpi so prepare_molecules takes the per-container branch."""
            return "mpi"

        def run(
            self,
            job_input: tuple[tuple[list[MolContainer], dict[str, object]], ...],
            func: Callable[..., None],
        ) -> None:
            """Capture the grouped jobs instead of dispatching them."""
            captured.extend(job_input)

        def end(self, job_manager: str) -> None:
            """No mpi universe to tear down."""

    monkeypatch.setattr(start, "Parallelizer", _StubParallelizer)

    start.prepare_molecules(
        {
            "source": str(src),
            "output_folder": str(tmp_path),
            "job_manager": "serial",
        }
    )

    assert len(captured) == 2
    for job_contnrs, _job_params in captured:
        contnr = job_contnrs[0]
        assert contnr.contnr_idx == 0
        assert contnr.mol_orig_frm_inp_smi.contnr_idx == 0

        # The manifestation: with a stale index on the pristine mol, the
        # container comes out of the filter empty.
        desalt_orig_smi([contnr], 1, "serial", None)
        durrant_lab_filters([contnr], 1, "serial", None)
        assert len(contnr.mols) == 1


def test_finalize_params_rejects_missing_source_file(tmp_path) -> None:
    # Regression: the existence check sat in an except clause around
    # os.path.abspath, which does not touch the filesystem and so never
    # raises. A typo'd filename fell through to load_smiles_file and surfaced
    # as a bare FileNotFoundError traceback.
    with pytest.raises(Exception, match="not found"):
        start.finalize_params(
            {
                "source": str(tmp_path / "typo.smi"),
                "output_folder": str(tmp_path),
                "add_pdb_output": False,
                "separate_output_files": False,
                "job_manager": "serial",
            }
        )


def test_prepare_molecules_rejects_missing_source_file(tmp_path) -> None:
    with pytest.raises(Exception, match="not found"):
        start.prepare_molecules(
            {"source": str(tmp_path / "nope.smi"), "output_folder": str(tmp_path)}
        )


def test_prepare_molecules_rejects_unsupported_source_extension(tmp_path) -> None:
    # Regression: an unrecognized extension used to put the source string
    # itself into smiles_data, which then failed while being unpacked as a
    # (smiles, name, properties) tuple.
    src = tmp_path / "molecules.txt"
    src.write_text("CCO\tethanol\n")
    with pytest.raises(Exception, match="extension"):
        start.prepare_molecules(
            {
                "source": str(src),
                "output_folder": str(tmp_path),
                "job_manager": "serial",
            }
        )


def test_prepare_molecules_rejects_unknown_json_parameter(tmp_path) -> None:
    src = tmp_path / "input.smi"
    src.write_text("CCO\tethanol\n")
    params_path = tmp_path / "params.json"
    params_path.write_text(json.dumps({"source": str(src), "bogus_flag": True}))
    with pytest.raises(Exception, match="Unrecognized parameter"):
        start.prepare_molecules({"json": str(params_path)})
