"""Unit tests for MyMol and MyConformer."""

import random

import numpy
import pytest
from rdkit import Chem

from gypsum_dl import MyMol


def test_mymol_from_rdkit_mol_sets_canonical_smiles() -> None:
    mol = MyMol.MyMol(Chem.MolFromSmiles("OCC"), "ethanol")
    # Asked through the accessor rather than read off can_smi: __init__ drops
    # the cache once the molecule has been sanitized, since sanitization can
    # replace the molecule the SMILES was taken from.
    assert mol.smiles() == "CCO"
    assert mol.name == "ethanol"


def test_mymol_smiles_describes_the_sanitized_molecule() -> None:
    # Regression: __init__ canonicalized the starter before
    # make_mol_frm_smiles_sanitze ran, and check_sanitization hands back a
    # modified copy when it applies the four-bond nitrogen fix. The cached
    # SMILES then described a molecule the object no longer held, and smiles()
    # kept returning it, so molecular identity (hashing, deduplication) and the
    # SMILES written to the SDF disagreed about the molecule.
    starter = Chem.MolFromSmiles("C[N](C)(C)C", sanitize=False)
    # Enough bookkeeping for the SMILES writer to run, but not a sanitization:
    # the neutral quaternary nitrogen has to reach MyMol uncorrected.
    starter.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(starter)
    assert "+" not in Chem.MolToSmiles(starter)

    mol = MyMol.MyMol(starter)

    assert mol.rdkit_mol is not None
    assert "+" in mol.smiles()
    assert mol.smiles() == Chem.MolToSmiles(
        mol.rdkit_mol, isomericSmiles=True, canonical=True
    )


def test_mymol_from_rdkit_mol_survives_smiles_conversion_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: when MolToSmiles raises for an RDKit-mol starter, __init__
    # used to reference the unbound local `smiles` and die with
    # UnboundLocalError instead of recording the failure.
    def boom(*args, **kwargs):
        raise ValueError("cannot canonicalize")

    monkeypatch.setattr(MyMol.Chem, "MolToSmiles", boom)
    mol = MyMol.MyMol(Chem.MolFromSmiles("OCC"), "ethanol")
    assert mol.can_smi is None
    assert mol.orig_smi == ""


def test_smiles_reports_failure_as_none_only(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: __init__ marked a failed canonicalization with False while
    # smiles() marked it with None, so the accessor had three possible return
    # types (str, None, False) and every consumer had to know which failure it
    # was looking at. There is now one non-string case.
    def boom(*args, **kwargs):
        raise ValueError("cannot canonicalize")

    monkeypatch.setattr(MyMol.Chem, "MolToSmiles", boom)

    from_rdkit_mol = MyMol.MyMol(Chem.MolFromSmiles("OCC"), "ethanol")
    assert from_rdkit_mol.smiles() is None

    from_smiles = MyMol.MyMol("CCO")
    assert from_smiles.smiles() is None


def test_molecules_with_unknown_smiles_are_not_interchangeable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: __hash__ hashed whatever smiles() handed back, so every
    # molecule whose canonical SMILES could not be computed shared one hash
    # and compared equal. The hash-based deduplication passes then kept one
    # such molecule and discarded the rest as copies of it.
    def boom(*args, **kwargs):
        raise ValueError("cannot canonicalize")

    monkeypatch.setattr(MyMol.Chem, "MolToSmiles", boom)
    first = MyMol.MyMol("CCO")
    second = MyMol.MyMol("CCCCCC")

    assert first.smiles() is None
    assert second.smiles() is None
    assert first != second
    assert hash(first) != hash(second)
    assert first == first
    assert len(dict.fromkeys([first, second])) == 2


def test_equality_against_a_non_molecule_is_false() -> None:
    # Regression: __eq__ asked the other operand for its hash, so comparing a
    # molecule against a bool sentinel compared hashes with hash(True) == 1.
    # A molecule whose canonical SMILES happened to hash to 1 would have
    # compared equal to True.
    mol = MyMol.MyMol("CCO")

    assert (mol == True) is False  # noqa: E712
    assert (mol != True) is True  # noqa: E712
    assert (mol == "CCO") is False
    assert (mol == None) is False  # noqa: E711


def test_smiles_noh_strips_explicit_hydrogens() -> None:
    assert MyMol.MyMol("CCO").smiles(True) == "CCO"


def test_smiles_noh_returns_none_when_there_is_no_molecule() -> None:
    # Regression: the noh == True branch had no failure guard. When rdkit_mol is
    # None (a molecule that failed sanitization or 3D optimization), the copy
    # and the deprotonation both give None and MolToSmiles raised a Boost
    # ArgumentError. set_all_rdkit_mol_props reaches this accessor from the main
    # process, after every expensive step, so that exception ended the run
    # before anything was written to disk.
    mol = MyMol.MyMol("CCO")
    mol.rdkit_mol = None

    assert mol.smiles(True) is None
    mol.set_all_rdkit_mol_props()  # must not raise


def test_smiles_is_cached() -> None:
    mol = MyMol.MyMol("CCO")
    assert mol.smiles() == "CCO"
    assert mol.smiles() == "CCO"
    assert mol.smiles(True) == mol.smiles(True)


def test_comparison_operators_follow_canonical_smiles() -> None:
    # Regression: every operator delegated to __hash__. str hashing is salted
    # per interpreter, so any sort that fell through to comparing molecules
    # (sorting (energy, MyMol) pairs whose energies tie) ordered them
    # differently on every run. Ordering follows the canonical SMILES itself
    # now, which is stable.
    a = MyMol.MyMol("CCO")
    b = MyMol.MyMol("OCC")
    c = MyMol.MyMol("CCC")
    assert a == b
    assert a != c
    assert a <= b
    assert a >= b
    assert (a < c) == (a.smiles() < c.smiles())
    assert (a > c) == (a.smiles() > c.smiles())
    assert sorted([a, c]) == [c, a]
    assert a.sort_key() == (0, "CCO")
    assert a is not None


def test_hash_collision_does_not_make_molecules_equal() -> None:
    # Regression: __eq__ compared hashes, so two distinct molecules whose
    # canonical SMILES happened to hash to the same value compared equal, and
    # the deduplication passes discarded one of them as a copy of the other.
    class CollidingMol(MyMol.MyMol):
        """A molecule that hashes the same as every other one of its kind.

        Collisions between real SMILES hashes are too rare to provoke from a
        test, so the collision is forced here instead.
        """

        def __hash__(self) -> int:
            """Collide with every other instance of this class."""
            return 1234

    first = CollidingMol("CCO")
    second = CollidingMol("CCC")

    assert hash(first) == hash(second)
    assert first != second
    assert (first == second) is False
    assert len(dict.fromkeys([first, second])) == 2


def test_ordering_ties_molecules_that_cannot_be_canonicalized(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Two molecules that both failed to canonicalize have not been shown to be
    # the same molecule, so neither is less than, greater than, or equal to the
    # other, and a stable sort leaves them in the order they arrived in. Under
    # the hash-based operators their order came from their memory addresses.
    def boom(*args, **kwargs):
        raise ValueError("cannot canonicalize")

    monkeypatch.setattr(MyMol.Chem, "MolToSmiles", boom)
    first = MyMol.MyMol("CCO")
    second = MyMol.MyMol("CCCCCC")

    assert first.sort_key() == second.sort_key()
    assert (first < second) is False
    assert (second < first) is False
    assert (first <= second) is False
    assert (first >= second) is False
    assert first != second
    assert sorted([first, second]) == [first, second]


def test_ordering_against_a_non_molecule_is_not_implemented() -> None:
    # Regression: the operators asked the other operand for its hash, so
    # comparing a molecule against an int compared canonical-SMILES hashes
    # with hash(5) == 5 and answered as though the comparison meant something.
    mol = MyMol.MyMol("CCO")

    assert mol.__lt__(5) is NotImplemented
    with pytest.raises(TypeError):
        mol < 5  # noqa: B015


def test_standardize_smiles_is_cached() -> None:
    mol = MyMol.MyMol("CCO")
    first = mol.standardize_smiles()
    assert first == mol.standardize_smiles()


def test_count_hyd_bnd_to_carb() -> None:
    assert MyMol.MyMol("CCO").count_hyd_bnd_to_carb() == 5


def test_get_idxs_of_nonaro_rng_atms_is_cached() -> None:
    mol = MyMol.MyMol("C1CCCCC1")
    rings = mol.get_idxs_of_nonaro_rng_atms()
    assert len(rings) == 1
    assert mol.get_idxs_of_nonaro_rng_atms() is rings


def test_get_idxs_of_nonaro_rng_atms_ignores_aromatic_rings() -> None:
    assert MyMol.MyMol("c1ccccc1").get_idxs_of_nonaro_rng_atms() == []


def test_chiral_center_helpers_are_cached() -> None:
    mol = MyMol.MyMol("CC(N)C(=O)O")
    unassigned = mol.chiral_cntrs_w_unasignd()
    assert len(unassigned) == 1
    assert mol.chiral_cntrs_w_unasignd() is unassigned
    assigned = mol.chiral_cntrs_only_asignd()
    assert assigned == []
    assert mol.chiral_cntrs_only_asignd() is assigned


def test_get_double_bonds_without_stereochemistry_finds_unspecified_bond() -> None:
    assert len(MyMol.MyMol("CC=CC").get_double_bonds_without_stereochemistry()) == 1


def test_get_double_bonds_without_stereochemistry_ignores_specified_bond() -> None:
    mol = MyMol.MyMol(Chem.MolFromSmiles(r"C/C=C/C"))
    assert mol.get_double_bonds_without_stereochemistry() == []


def test_remove_bizarre_substruc_flags_carbanion() -> None:
    mol = MyMol.MyMol("CC[CH2-]")
    assert mol.remove_bizarre_substruc() is True
    assert mol.remove_bizarre_substruc() is True


def test_remove_bizarre_substruc_allows_normal_molecule() -> None:
    mol = MyMol.MyMol("CCO")
    assert mol.remove_bizarre_substruc() is False
    assert mol.remove_bizarre_substruc() is False


def test_remove_bizarre_substruc_survives_non_string_can_smi() -> None:
    # Regression (M8): can_smi is not a string once canonicalization has
    # failed (None now, False in older versions); `s in self.can_smi` then
    # raised TypeError ("argument of type 'bool'/'NoneType' is not iterable"),
    # which became a hang under multiprocessing. The method must return a bool
    # instead, whichever non-string marker it finds.
    mol = MyMol.MyMol("CCO")
    mol.can_smi = False
    result = mol.remove_bizarre_substruc()
    assert isinstance(result, bool)

    mol2 = MyMol.MyMol("CCO")
    mol2.can_smi = None
    assert isinstance(mol2.remove_bizarre_substruc(), bool)


def test_get_frags_of_orig_smi_single_fragment_returns_self() -> None:
    mol = MyMol.MyMol("CCO")
    assert mol.get_frags_of_orig_smi() == [mol]


def test_get_frags_of_orig_smi_splits_salts() -> None:
    assert len(MyMol.MyMol("CCO.CC").get_frags_of_orig_smi()) == 2


def test_inherit_contnr_props() -> None:
    source = MyMol.MyMol("CCO", "ethanol")
    source.contnr_idx = 3
    target = MyMol.MyMol("CCC")
    target.inherit_contnr_props(source)
    assert target.contnr_idx == 3
    assert target.name == "ethanol"
    assert target.orig_smi == source.orig_smi


def test_set_all_rdkit_mol_props_records_genealogy() -> None:
    mol = MyMol.MyMol("CCO", "ethanol")
    mol.mol_props["activity"] = 1.5
    mol.genealogy = ["CCO (source)"]
    mol.set_all_rdkit_mol_props()
    assert mol.rdkit_mol.GetProp("SMILES") == "CCO"
    assert mol.rdkit_mol.GetProp("activity") == "1.5"
    assert mol.rdkit_mol.GetProp("Genealogy") == "CCO (source)"
    assert mol.rdkit_mol.GetProp("_Name") == "ethanol"


def test_set_all_rdkit_mol_props_omits_an_unknown_smiles() -> None:
    # Regression (bug 5): when smiles(True) failed it returned None, which
    # str()'d to the literal "None" and was written to the SDF as the molecule's
    # SMILES. A failed calculation now leaves the field absent instead.
    mol = MyMol.MyMol("CCO", "ethanol")
    mol.can_smi_noh = None
    mol.genealogy = ["CCO (source)"]

    mol.set_all_rdkit_mol_props()

    assert not mol.rdkit_mol.HasProp("SMILES")
    # The remaining properties are still written.
    assert mol.rdkit_mol.GetProp("Genealogy") == "CCO (source)"
    assert mol.rdkit_mol.GetProp("_Name") == "ethanol"


def test_set_rdkit_mol_prop_writes_once_and_tolerates_none() -> None:
    # Regression (B14): set_rdkit_mol_prop used to SetProp three times, the
    # first two unguarded against rdkit_mol being None. Collapsed to a single
    # guarded write: the property is set, and a None rdkit_mol no longer raises.
    mol = MyMol.MyMol("CCO")
    mol.set_rdkit_mol_prop("activity", 1.5)
    assert mol.rdkit_mol.GetProp("activity") == "1.5"

    mol.rdkit_mol = None
    mol.set_rdkit_mol_prop("activity", 1.5)  # must not raise


def test_myconformer_moltomolblock_logs_block(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Regression (B11): MolToMolBlock referenced self.mol_copy (nonexistent) and
    # self.conformer (a method, not an attribute), so it raised instead of
    # logging. It must now emit a molblock.
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    mol.conformers[0].MolToMolBlock()
    assert "RDKit" in capsys.readouterr().out


def test_make_first_3d_conf_no_min_is_idempotent() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    assert len(mol.conformers) == 1
    mol.make_first_3d_conf_no_min()
    assert len(mol.conformers) == 1


def test_make_first_3d_conf_no_min_survives_failed_reprotanation(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression (M13): try_reprotanation can return None. That None was
    # assigned straight to rdkit_mol, and add_conformers -> MyConformer then did
    # copy.deepcopy(None).RemoveAllConformers(), raising AttributeError inside
    # pick_lowest_enrgy_mols. The method must bail out quietly instead.
    monkeypatch.setattr(MyMol.MOH, "try_reprotanation", lambda mol: None)
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    assert mol.conformers == []


def test_myconformer_with_none_rdkit_mol_is_marked_failed() -> None:
    # Regression (M13): constructing a MyConformer from a MyMol whose rdkit_mol
    # is None must not crash (deepcopy(None).RemoveAllConformers()). It should
    # instead flag itself failed (mol is False) so add_conformers skips it.
    mol = MyMol.MyMol("CCO")
    mol.rdkit_mol = None
    conf = MyMol.MyConformer(mol)
    assert conf.mol is False


def test_add_conformers_sorts_by_energy() -> None:
    mol = MyMol.MyMol("CCCCCC")
    # `MyConformer.rmsd_to_me` rebuilds the molecule from SMILES and
    # reprotonates it, so it only matches conformers whose parent already
    # carries explicit hydrogens. The pipeline guarantees this by calling
    # `make_first_3d_conf_no_min` before any RMSD-based pruning; do the same
    # here rather than embedding an implicit-hydrogen molecule.
    mol.make_first_3d_conf_no_min()
    mol.add_conformers(3, 0.1, True)
    energies = [conf.energy for conf in mol.conformers]
    assert len(energies) >= 1
    assert energies == sorted(energies)


def test_load_conformers_into_rdkit_mol() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    mol.load_conformers_into_rdkit_mol()
    assert mol.rdkit_mol.GetNumConformers() == 1


def test_load_conformers_into_rdkit_mol_assigns_distinct_ids() -> None:
    # Regression: AddConformer was called without assignId, so each conformer
    # kept its source id. Every freshly embedded MyConformer holds exactly one
    # conformer with id 0, so a multi-conformer molecule ended up with several
    # conformers sharing id 0. SDWriter and MolToPDBFile resolve confId=-1 to
    # the first match, so the same coordinates would be written repeatedly.
    mol = MyMol.MyMol("CCCCCC")
    mol.make_first_3d_conf_no_min()
    mol.conformers.append(MyMol.MyConformer(mol))

    mol.load_conformers_into_rdkit_mol()

    ids = sorted(conf.GetId() for conf in mol.rdkit_mol.GetConformers())
    assert ids == [0, 1]


def test_conformer_coords_and_energy() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    assert conf.coords().shape[0] == conf.mol.GetNumAtoms()
    assert isinstance(conf.get_energy(), float)


def test_conformer_minimize_is_idempotent() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    conf.minimize()
    energy = conf.energy
    conf.minimize()
    assert conf.energy == energy


def test_conformer_write_pdb_file(tmp_path) -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    out = tmp_path / "conf.pdb"
    mol.conformers[0].write_pdb_file(str(out))
    # Small molecules carry no residue information, so RDKit emits HETATM.
    assert "HETATM" in out.read_text()


def test_conformer_rmsd_between_identical_conformers_is_zero() -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    duplicate = MyMol.MyConformer(mol, conf.conformer())
    assert conf.rmsd_to_me(duplicate) == pytest.approx(0.0, abs=1e-6)


def test_conformer_rmsd_is_infinite_when_the_smiles_is_unavailable() -> None:
    # Regression: rmsd_to_me handed self.smiles straight to MolFromSmiles, but
    # MyMol.smiles() reports failure as None, and MolFromSmiles(None) raises.
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    duplicate = MyMol.MyConformer(mol, conf.conformer())
    conf.smiles = None

    assert conf.rmsd_to_me(duplicate) == float("inf")


def test_eliminate_similar_conformers_keeps_both_when_rmsd_cannot_be_computed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: rmsd_to_me chained check_sanitization and try_reprotanation
    # and then called AddConformer on the result without a None check, unlike
    # every other consumer of those helpers. The AttributeError escaped
    # add_conformers and the ring-conformer and minimization steps, reaching
    # parallelizer.run_one, which printed the traceback and dropped the whole
    # molecule. Keeping both conformers is the conservative outcome: nothing is
    # deduplicated on an RMSD that was never computed.
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    mol.conformers.append(MyMol.MyConformer(mol, mol.conformers[0].conformer()))
    monkeypatch.setattr(MyMol.MOH, "try_reprotanation", lambda amol: None)

    mol.eliminate_structurally_similar_conformers(0.1)

    assert len(mol.conformers) == 2


def test_eliminate_similar_conformers_still_drops_duplicates() -> None:
    # Control for the test above: two copies of the same conformer do collapse
    # when the RMSD is computable, so the infinity fallback is not masking the
    # deduplication entirely.
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    mol.conformers.append(MyMol.MyConformer(mol, mol.conformers[0].conformer()))

    mol.eliminate_structurally_similar_conformers(0.1)

    assert len(mol.conformers) == 1


def test_coord_3d_err_warning_is_logged(capsys: pytest.CaptureFixture[str]) -> None:
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    mol.conformers[0].coord_3d_err_warning(None)
    assert "WARNING" in capsys.readouterr().out


def _first_conformer_positions(smiles: str, seed: int) -> numpy.ndarray:
    """Embed a molecule's first conformer under a fixed Python seed.

    Written as a helper so the two seeded embeddings the test compares are
    built the same way, each starting from a freshly seeded Python generator.

    Args:
        smiles: SMILES string of the molecule to embed.
        seed: Seed handed to the random module before embedding.

    Returns:
        The coordinates of the resulting first conformer.
    """
    random.seed(seed)
    mol = MyMol.MyMol(smiles)
    mol.make_first_3d_conf_no_min()
    return mol.conformers[0].coords()


def test_conformer_coordinates_are_reproducible_under_a_fixed_seed() -> None:
    # Regression: MyConformer never set params.randomSeed, and RDKit's
    # embedding keeps a generator of its own that random.seed() cannot reach.
    # Its default (-1) picks a fresh seed per call, so --random_seed left the
    # output coordinates varying between otherwise identical serial runs. The
    # chain is long enough to have several accessible geometries.
    first = _first_conformer_positions("CCCCCCCCO", 202609)
    second = _first_conformer_positions("CCCCCCCCO", 202609)

    assert numpy.allclose(first, second)


def test_second_embed_fallback_passes_a_seed(monkeypatch: pytest.MonkeyPatch) -> None:
    # Regression: the legacy rescue embedder takes no EmbedParameters, so it
    # kept RDKit's unseeded default even once the primary call was seeded,
    # leaving any molecule rescued by --second_embed nonreproducible.
    calls: list[dict[str, object]] = []

    def recording_embed(mol: Chem.Mol, *args: object, **kwargs: object) -> int:
        """Record each embedding attempt and report that it produced nothing.

        Returning failure without adding a conformer walks MyConformer through
        both retries, so the legacy fallback is the last call recorded.

        Args:
            mol: The molecule handed to RDKit.
            *args: Positional arguments, which carry the EmbedParameters.
            **kwargs: Keyword arguments, which carry the legacy call's seed.

        Returns:
            RDKit's failure code.
        """
        calls.append(dict(kwargs))
        return -1

    monkeypatch.setattr(MyMol.AllChem, "EmbedMolecule", recording_embed)

    mol = MyMol.MyMol("CCO")
    random.seed(202609)
    MyMol.MyConformer(mol, second_embed=True)

    assert calls, "EmbedMolecule was never called"
    seed = calls[-1].get("randomSeed")
    assert isinstance(seed, int)
    assert seed > 0


def test_standardize_smiles_survives_an_unknown_noh_smiles(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    # Regression: the error handler built its message from smiles(True), which
    # is None once deprotonation has failed, so a molvs failure was replaced by
    # a TypeError raised inside the handler itself. That aborted the run after
    # the SDF had already been written, during PDB output.
    mol = MyMol.MyMol("CCO")
    mol.can_smi_noh = None

    def boom(smiles: str) -> str:
        raise ValueError("cannot standardize")

    monkeypatch.setattr(MyMol, "ssmiles", boom)

    assert mol.standardize_smiles() == mol.smiles()
    assert "Could not standardize" in capsys.readouterr().out


def test_conformer_energy_is_infinite_when_the_force_field_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: a failed force field recorded the magic number 9999, which a
    # strained or large ligand can legitimately exceed, so an unscored
    # conformer could outrank a real one everywhere energies are sorted.
    mol = MyMol.MyMol("CCO")

    def boom(mol_to_score: object) -> object:
        raise ValueError("no UFF parameters")

    monkeypatch.setattr(MyMol.AllChem, "UFFGetMoleculeForceField", boom)

    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    assert conf.energy == float("inf")

    conf.minimized = False
    conf.minimize()
    assert conf.energy == float("inf")
