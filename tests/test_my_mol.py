"""Unit tests for MyMol and MyConformer."""

import copy
import pickle
import random

import numpy
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

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


def test_mymol_reports_no_smiles_when_sanitization_fails(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: __init__ canonicalized an RDKit-mol starter before
    # sanitizing it and cleared that cache only when sanitization succeeded. A
    # molecule check_sanitization rejected was therefore left with rdkit_mol
    # None but a healthy-looking cached SMILES, which defeats every caller
    # that treats a non-string SMILES as "molecule of unknown identity":
    # uniq_mols_in_list deduplicated it against unrelated molecules, and
    # remove_highly_charged_molecules would hand None to
    # Chem.GetFormalCharge.
    monkeypatch.setattr(MyMol.MOH, "check_sanitization", lambda mol: None)

    from_rdkit_mol = MyMol.MyMol(Chem.MolFromSmiles("CCO"), "ethanol")

    assert from_rdkit_mol.rdkit_mol is None
    assert from_rdkit_mol.can_smi is None
    assert from_rdkit_mol.smiles() is None

    # The SMILES-starter path never had a pre-sanitization SMILES to keep, so
    # it still discovers the failure in smiles() itself.
    from_smiles = MyMol.MyMol("CCO", "ethanol")

    assert from_smiles.rdkit_mol is None
    assert from_smiles.smiles() is None


def test_mymol_does_not_alias_the_callers_rdkit_mol() -> None:
    # Regression: __init__ stored the starter itself, and check_sanitization
    # hands the same object back whenever the molecule sanitizes on the first
    # pass, so a MyMol and the code that built it shared one mutable molecule.
    # load_conformers_into_rdkit_mol clears the conformer set and, with no
    # conformers to load, computes 2D coordinates, which would appear on the
    # caller's molecule underneath it.
    starter = Chem.MolFromSmiles("CCO")
    assert starter.GetNumConformers() == 0

    mol = MyMol.MyMol(starter)
    assert mol.rdkit_mol is not starter

    mol.load_conformers_into_rdkit_mol()

    assert mol.rdkit_mol.GetNumConformers() == 1
    assert starter.GetNumConformers() == 0


def test_mymol_conformers_do_not_overwrite_the_callers_coordinates() -> None:
    # The same aliasing seen from the other direction: a caller that keeps its
    # own molecule, coordinates and all, must not have them replaced by the 3D
    # conformer this object generates.
    starter = Chem.AddHs(Chem.MolFromSmiles("CCO"))
    AllChem.Compute2DCoords(starter)
    coords_before = starter.GetConformer().GetPositions().tolist()

    mol = MyMol.MyMol(starter)
    mol.add_conformers(1, 1e60, False)
    mol.load_conformers_into_rdkit_mol()

    assert starter.GetNumConformers() == 1
    assert starter.GetConformer().GetPositions().tolist() == coords_before


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


def test_remove_bizarre_substruc_consults_the_structure_not_the_smiles() -> None:
    # Regression: a substring pass over orig_smi, orig_smi_deslt, and can_smi
    # ran ahead of the substructure matching, so a prohibited pattern that
    # merely appeared in one of those strings rejected the molecule even
    # though the molecule itself did not contain it. Only the structure should
    # decide.
    mol = MyMol.MyMol("CCO")
    assert mol.rdkit_mol is not None
    mol.orig_smi = "[C-]#[O+].CCO"
    mol.orig_smi_deslt = "[C-]#[O+].CCO"

    assert mol.remove_bizarre_substruc() is False


@pytest.mark.parametrize(
    "smiles",
    [
        "CC(=C)O",  # C(=[CH2])[OH], a terminal enol.
        "C=C(O)O",  # C=C([OH])[OH], a geminal vinyl diol.
        "[C-]#[O+]",  # [C-], a carbanion.
    ],
)
def test_remove_bizarre_substruc_matches_patterns_the_substring_pass_missed(
    smiles: str,
) -> None:
    # These patterns are SMARTS and cannot appear verbatim in the SMILES
    # strings the removed substring pass compared against, so the substructure
    # matching is the only thing enforcing them. Assert the molecules sanitize
    # first, since an unbuildable molecule is reported as bizarre for an
    # unrelated reason.
    mol = MyMol.MyMol(smiles)
    assert mol.rdkit_mol is not None
    assert mol.remove_bizarre_substruc() is True


def test_get_frags_of_orig_smi_single_fragment_returns_the_wrapped_mol() -> None:
    # The single-fragment branch used to return [self], so the element type of
    # the list depended on the fragment count. The desalter reads
    # GetNumHeavyAtoms off these elements on its multi-fragment path, and was
    # safe only because the single-fragment case returns before reaching it.
    mol = MyMol.MyMol("CCO")
    frags = mol.get_frags_of_orig_smi()
    assert frags == [mol.rdkit_mol]
    assert frags[0].GetNumHeavyAtoms() == 3


def test_get_frags_of_orig_smi_splits_salts() -> None:
    frags = MyMol.MyMol("CCO.CC").get_frags_of_orig_smi()
    assert len(frags) == 2
    assert all(hasattr(frag, "GetNumHeavyAtoms") for frag in frags)


def test_get_frags_of_orig_smi_caches_the_list_it_returns() -> None:
    mol = MyMol.MyMol("CCO")
    assert mol.get_frags_of_orig_smi() is mol.get_frags_of_orig_smi()


def test_inherit_contnr_props() -> None:
    source = MyMol.MyMol("CCO", "ethanol")
    source.contnr_idx = 3
    target = MyMol.MyMol("CCC")
    target.inherit_contnr_props(source)
    assert target.contnr_idx == 3
    assert target.name == "ethanol"
    assert target.orig_smi == source.orig_smi


def test_inherit_contnr_props_accepts_an_extracted_mapping() -> None:
    # The tautomer step ships these fields to its workers in place of the
    # container, so the same field list has to apply on the far side.
    target = MyMol.MyMol("CCC")
    target.inherit_contnr_props(
        {
            "contnr_idx": 4,
            "name": "ethanol",
            "orig_smi": "CCO",
            "orig_smi_deslt": "CCO",
            "orig_smi_canonical": "CCO",
        }
    )
    assert target.contnr_idx == 4
    assert target.name == "ethanol"
    assert target.orig_smi == "CCO"
    assert target.orig_smi_canonical == "CCO"


def test_orig_smi_canonical_exists_on_every_molecule() -> None:
    # It used to be set only by MolContainer.add_smiles, so which fields a
    # MyMol carried depended on which code path built it, and a variant built
    # by inherit_contnr_props (every tautomer, enantiomer, and cis/trans form)
    # never acquired it at all.
    assert MyMol.MyMol("CCO").orig_smi_canonical is None


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


def test_set_all_rdkit_mol_props_prefers_the_computed_smiles() -> None:
    # Regression: SMILES was written before mol_props was iterated, and
    # mol_props carries the input file's tags, so an input SDF tag named SMILES
    # replaced the prepared variant's own SMILES. Every variant of that input
    # then advertised the same (possibly salted, possibly mis-protonated) input
    # string, which collapses the variants for anyone keying off the field.
    mol = MyMol.MyMol("CCO", "ethanol")
    mol.mol_props["SMILES"] = "CCO.[Na+]"

    mol.set_all_rdkit_mol_props()

    assert mol.rdkit_mol.GetProp("SMILES") == "CCO"


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


def test_myconformer_has_an_energy_even_when_it_fails() -> None:
    # Regression (F6): the failure paths returned (or fell past) the block that
    # assigns self.energy, so the attribute never existed. Two callers read it
    # without checking .mol first (the minimization and ring-conformer steps),
    # and add_conformers sorts on it, so a failed conformer raised
    # AttributeError inside a worker and the molecule was dropped with no
    # explanation. Infinity also keeps a scoreless conformer from outranking a
    # real one.
    mol = MyMol.MyMol("CCO")
    mol.rdkit_mol = None
    conf = MyMol.MyConformer(mol)
    assert conf.energy == float("inf")


def test_failed_myconformer_is_fully_initialized() -> None:
    # Regression: the rdkit_mol-is-None path returned before self.minimized
    # and self.ids_hvy_atms were assigned, so minimize() and align_to_me()
    # raised AttributeError on such an object. minimize() was worse than that:
    # it handed False to the force field, and its own handler then raised
    # again on Chem.MolToSmiles(False), so the exception escaped the except
    # block meant to absorb it. Minimize3D.parallel_minit builds exactly this
    # object and stores it as the molecule's only conformer.
    mol = MyMol.MyMol("CCO")
    mol.rdkit_mol = None
    conf = MyMol.MyConformer(mol)

    assert conf.mol is False
    assert conf.minimized is True
    assert conf.ids_hvy_atms == []
    conf.minimize()
    assert conf.energy == float("inf")


def test_add_conformers_sorts_by_energy() -> None:
    mol = MyMol.MyMol("CCCCCC")
    # The pipeline always calls `make_first_3d_conf_no_min` before any
    # RMSD-based pruning, which is what puts explicit hydrogens on the parent;
    # do the same here rather than embedding an implicit-hydrogen molecule.
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


def test_conformer_rmsd_is_infinite_when_the_molecule_is_unavailable() -> None:
    # A conformer that failed to embed carries mol is False, so there is no
    # graph to compare coordinates in. Infinity is the conservative answer:
    # eliminate_structurally_similar_conformers deduplicates on
    # rmsd <= cutoff, so nothing is discarded on an uncomputable RMSD.
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    duplicate = MyMol.MyConformer(mol, conf.conformer())
    conf.mol = False

    assert conf.rmsd_to_me(duplicate) == float("inf")
    assert duplicate.rmsd_to_me(conf) == float("inf")


def _hydrogens_first(conf: "MyMol.MyConformer") -> None:
    """Renumber a conformer's molecule so its hydrogens come first.

    rmsd_to_me used to rebuild the molecule from the canonical SMILES, whose
    atom order puts heavy atoms first, and then index the conformer's
    coordinates by that rebuilt order. Any molecule whose own order already
    happens to be heavy-atoms-first hides the mismatch, so the test has to
    supply one that is not (an SDF input, for instance, can carry hydrogens
    anywhere in the block).

    Args:
        conf: The conformer to renumber, modified in place along with the
            heavy-atom index list align_to_me reads.
    """
    order = sorted(
        range(conf.mol.GetNumAtoms()),
        key=lambda idx: conf.mol.GetAtomWithIdx(idx).GetAtomicNum(),
    )
    conf.mol = Chem.RenumberAtoms(conf.mol, order)
    conf.ids_hvy_atms = [
        a.GetIdx() for a in conf.mol.GetAtoms() if a.GetAtomicNum() != 1
    ]


def test_conformer_rmsd_uses_the_coordinates_own_atom_ordering() -> None:
    # Regression: rmsd_to_me built its comparison molecule from
    # self.smiles (the canonical SMILES), whose atom order is unrelated to
    # that of the molecule the coordinates came from. Deprotonating that
    # rebuilt molecule then dropped whichever atoms were hydrogens in the
    # rebuilt ordering, so the reported number was an RMSD over a mixed subset
    # of heavy atoms and hydrogens instead of a heavy-atom RMSD. Here the
    # heavy atoms are displaced and the hydrogens are not, so the two answers
    # are far apart: 5 A over the heavy atoms, or 0 A over the hydrogens the
    # old code compared in their place.
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    conf = mol.conformers[0]
    _hydrogens_first(conf)

    shifted = Chem.Conformer(conf.conformer())
    for atom in conf.mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            continue
        pos = shifted.GetAtomPosition(atom.GetIdx())
        shifted.SetAtomPosition(atom.GetIdx(), Point3D(pos.x + 5.0, pos.y, pos.z))

    other = copy.deepcopy(conf)
    other.conformer(shifted)

    assert conf.rmsd_to_me(other) == pytest.approx(5.0, abs=1e-6)


def test_eliminate_similar_conformers_keeps_both_when_rmsd_cannot_be_computed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: rmsd_to_me called AddConformer and GetConformerRMS on the
    # result of a MolObjectHandling helper without a None check, unlike every
    # other consumer of those helpers. The AttributeError escaped
    # add_conformers and the ring-conformer and minimization steps, reaching
    # parallelizer.run_one, which printed the traceback and dropped the whole
    # molecule. Keeping both conformers is the conservative outcome: nothing is
    # deduplicated on an RMSD that was never computed.
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()
    mol.conformers.append(MyMol.MyConformer(mol, mol.conformers[0].conformer()))
    monkeypatch.setattr(MyMol.MOH, "try_deprotanation", lambda amol: None)

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


def test_unfilled_cache_marker_refuses_truthiness() -> None:
    # Regression: the "not computed yet" state used to be an empty string, so
    # the read sites had to distinguish it from a legitimately falsy cached
    # value (remove_bizarre_substruc caches False,
    # get_idxs_of_nonaro_rng_atms caches an empty list, smiles() caches None).
    # A truthiness test against the marker is always the wrong question, so it
    # must fail loudly rather than pick a branch.
    with pytest.raises(TypeError):
        bool(MyMol.UNSET)

    with pytest.raises(TypeError):
        if MyMol.UNSET:
            pass


def test_unfilled_cache_marker_survives_pickling_and_copying() -> None:
    # The parallelizer pickles MyMol objects out to workers and back at every
    # stage, and several steps deep-copy them. If the marker were rebuilt as a
    # separate instance, `is UNSET` would report every unfilled cache as
    # filled on the far side of the round trip, and the accessors would hand
    # back the marker itself as a result.
    assert pickle.loads(pickle.dumps(MyMol.UNSET)) is MyMol.UNSET
    assert copy.deepcopy(MyMol.UNSET) is MyMol.UNSET
    assert copy.copy(MyMol.UNSET) is MyMol.UNSET

    mol = MyMol.MyMol("CCO")
    revived = pickle.loads(pickle.dumps(mol))

    assert revived.can_smi is MyMol.UNSET
    assert revived.bizarre_substruct is MyMol.UNSET
    assert revived.smiles() == "CCO"


def test_cached_false_substructure_verdict_is_reused() -> None:
    # A cached False is a real answer, not an empty cache. This molecule does
    # contain a prohibited substructure, so a cache that was not consulted
    # would return True.
    mol = MyMol.MyMol("C=C(O)O")
    assert mol.remove_bizarre_substruc() is True

    mol.bizarre_substruct = False
    assert mol.remove_bizarre_substruc() is False


def test_cached_empty_ring_list_is_reused() -> None:
    # Same shape as the verdict cache: an empty list is the correct answer for
    # a molecule with no nonaromatic rings, so it has to be distinguishable
    # from a cache that has not been filled. Cyclohexane has one such ring.
    mol = MyMol.MyMol("C1CCCCC1")
    assert len(mol.get_idxs_of_nonaro_rng_atms()) == 1

    mol.nonaro_ring_atom_idx = []
    assert mol.get_idxs_of_nonaro_rng_atms() == []


def test_smiles_does_not_change_when_3d_coordinates_are_added() -> None:
    # Regression: make_first_3d_conf_no_min replaces rdkit_mol with the
    # AddHs version, and add_conformers then reaches MyConformer, which reads
    # smiles(). An unfilled cache was therefore filled from the
    # hydrogen-added molecule, so the canonical SMILES came out with explicit
    # hydrogens. Whether that happened depended only on whether something
    # earlier in the run had happened to call smiles() first.
    mol = MyMol.MyMol("CCO")
    mol.make_first_3d_conf_no_min()

    assert mol.conformers
    assert mol.smiles() == "CCO"
    assert "[H]" not in mol.smiles()


def test_smiles_is_the_same_whether_or_not_it_was_read_before_3d() -> None:
    # The two orderings have to agree: a variant that was deduplicated (which
    # calls smiles() through __hash__) before the 3D step must end up with the
    # same identity key, and the same reported SMILES, as one that was not.
    read_first = MyMol.MyMol("CC(=O)O")
    before = read_first.smiles()
    read_first.make_first_3d_conf_no_min()

    read_later = MyMol.MyMol("CC(=O)O")
    read_later.make_first_3d_conf_no_min()

    assert read_first.smiles() == before
    assert read_later.smiles() == before
    assert read_first.conformers[0].smiles == before
    assert read_later.conformers[0].smiles == before


def test_failed_reprotanation_leaves_the_smiles_cache_alone(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The cache is pinned only on the path that actually swaps rdkit_mol, so a
    # molecule that bails out early is left to canonicalize on demand as
    # before, and no canonicalization warning is emitted on its behalf.
    monkeypatch.setattr(MyMol.MOH, "try_reprotanation", lambda mol: None)

    mol = MyMol.MyMol("CCO")
    assert mol.can_smi is MyMol.UNSET

    mol.make_first_3d_conf_no_min()

    assert mol.conformers == []
    assert mol.can_smi is MyMol.UNSET


def test_index_caches_survive_the_hydrogen_swap() -> None:
    # AddHs appends the hydrogens, so the ring indices computed on the
    # hydrogen-free molecule still name the same atoms and are kept as is.
    mol = MyMol.MyMol("OC1CCCCC1")
    rings = mol.get_idxs_of_nonaro_rng_atms()
    mol.make_first_3d_conf_no_min()
    assert mol.nonaro_ring_atom_idx is rings


def test_index_caches_are_dropped_when_the_swap_renumbers_atoms(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: the ring and chiral-center caches survived the swap to the
    # hydrogen-added molecule unconditionally, so they were right only because
    # AddHs happens not to renumber. A renumbering swap left ring "members"
    # pointing at hydrogens.
    mol = MyMol.MyMol("OC1CCCCC1")
    mol.get_idxs_of_nonaro_rng_atms()
    mol.chiral_cntrs_w_unasignd()
    mol.chiral_cntrs_only_asignd()

    real_reprotanation = MyMol.MOH.try_reprotanation

    def renumbered(rdkit_mol: Chem.Mol) -> Chem.Mol:
        with_hs = real_reprotanation(rdkit_mol)
        order = list(reversed(range(with_hs.GetNumAtoms())))
        return Chem.RenumberAtoms(with_hs, order)

    monkeypatch.setattr(MyMol.MOH, "try_reprotanation", renumbered)
    mol.make_first_3d_conf_no_min()

    rings = mol.get_idxs_of_nonaro_rng_atms()
    assert len(rings) == 1
    assert [mol.rdkit_mol.GetAtomWithIdx(i).GetSymbol() for i in rings[0]] == ["C"] * 6
