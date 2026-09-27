"""Regression tests for SMILES desalting."""

import pytest
from rdkit import Chem

from gypsum_dl import MyMol
from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.steps.smiles import DeSaltOrigSmiles
from gypsum_dl.steps.smiles.DeSaltOrigSmiles import desalt_orig_smi, desalter


def test_desalter_does_not_alias_original_genealogy() -> None:
    # Regression (bug 7): desalter shared the pristine molecule's genealogy list
    # (`new_mol.genealogy = contnr.mol_orig_frm_inp_smi.genealogy`), so the
    # later append in desalt_orig_smi mutated the container's original molecule.
    contnr = MolContainer("CCCCO.[Na+]", "salt", 0, {})
    orig = contnr.mol_orig_frm_inp_smi.genealogy

    new_mol = desalter(contnr)

    # desalter appends the source entry to this very list, so snapshot it
    # afterward; what matters here is that the copy handed to the fragment is a
    # separate object.
    before = list(orig)
    assert new_mol.genealogy is not orig
    new_mol.genealogy.append("marker")
    assert contnr.mol_orig_frm_inp_smi.genealogy == before


def test_desalter_seeds_the_source_smiles_on_the_single_fragment_path() -> None:
    # Regression (bug 15): only the ionization-failure fallback ever wrote a
    # "(source)" entry, so a molecule that sailed through every step had a
    # genealogy starting at "(protonated)" and no record of its input SMILES.
    contnr = MolContainer("CCO", "ethanol", 0, {})

    mol = desalter(contnr)

    assert mol.genealogy[0] == "CCO (source)"


def test_desalter_breaks_a_size_tie_on_input_order() -> None:
    # Regression (bug 4): fragments were stored in a dict keyed by heavy-atom
    # count, so equal-sized fragments overwrote each other and the tie went to
    # whichever fragment came last. Benzene and hexane both have six heavy
    # atoms, so the kept fragment must now follow the order of the input.
    first = desalter(MolContainer("c1ccccc1.CCCCCC", "tie", 0, {}))
    assert Chem.CanonSmiles(first.smiles()) == Chem.CanonSmiles("c1ccccc1")

    reversed_input = desalter(MolContainer("CCCCCC.c1ccccc1", "tie_reversed", 1, {}))
    assert Chem.CanonSmiles(reversed_input.smiles()) == Chem.CanonSmiles("CCCCCC")


def test_desalter_seeds_the_source_smiles_on_the_desalted_path() -> None:
    # Regression (bug 15): the desalted branch copied the pristine molecule's
    # genealogy, which was empty, so the salted input SMILES was never
    # recorded either.
    contnr = MolContainer("CCCCO.[Na+]", "salt", 0, {})

    mol = desalter(contnr)

    assert mol.genealogy[0] == "CCCCO.[Na+] (source)"


def test_desalt_orig_smi_keeps_the_source_record_after_rebuilding_the_container() -> (
    None
):
    # Regression (bug 15): update_orig_smi replaces mol_orig_frm_inp_smi with a
    # fresh molecule whose genealogy is empty, and the ionization step seeds
    # every variant from that object. The source and desalt entries have to
    # survive the rebuild, or desalted molecules lose their provenance
    # downstream.
    salted = MolContainer("CCCCO.[Na+]", "salted", 0, {})

    desalt_orig_smi([salted], 1, "serial", None)

    genealogy = salted.mol_orig_frm_inp_smi.genealogy
    assert genealogy[0] == "CCCCO.[Na+] (source)"
    assert genealogy[-1].endswith("(desalted)")


def test_desalt_orig_smi_pairs_each_container_with_its_own_mol() -> None:
    # Regression (bug 6): pairing desalted mols to containers by a shared index
    # into a strip_none'd list could misalign them. Each container must receive
    # the fragment derived from its own input.
    salted = MolContainer("CCCCO.[Na+]", "salted", 0, {})
    clean = MolContainer("c1ccccc1", "benzene", 1, {})
    contnrs = [salted, clean]

    desalt_orig_smi(contnrs, 1, "serial", None)

    assert len(salted.mols) == 1
    assert len(clean.mols) == 1
    # Largest fragment of the salt is the butanol, and the salt was stripped.
    assert "Na" not in salted.mols[0].smiles()
    assert salted.orig_smi_deslt == salted.mols[0].orig_smi_deslt
    # The already-clean container is untouched apart from gaining its own mol.
    assert clean.mols[0].smiles() == clean.orig_smi_canonical


def test_desalt_orig_smi_isolates_a_container_whose_desalting_raises(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: desalting ran bare in the main process
    # (`[desalter(x) for x in contnrs]`), so a single molecule whose largest
    # fragment will not sanitize on its own aborted the whole run at the first
    # step, before any output existed.
    salted = MolContainer("CCCCO.[Na+]", "salted", 0, {})
    clean = MolContainer("c1ccccc1", "benzene", 1, {})

    real_desalter = DeSaltOrigSmiles.desalter

    def exploding_desalter(contnr: MolContainer) -> MyMol.MyMol:
        """Stand in for a fragment that cannot be sanitized in isolation.

        Args:
            contnr: The container being desalted.

        Returns:
            The desalted molecule, for every container but the salted one.
        """
        if contnr.name == "salted":
            raise ValueError("fragment will not sanitize on its own")
        return real_desalter(contnr)

    monkeypatch.setattr(DeSaltOrigSmiles, "desalter", exploding_desalter)

    desalt_orig_smi([salted, clean], 1, "serial", None)

    # The healthy container is unaffected.
    assert len(clean.mols) == 1
    assert clean.mols[0].smiles() == clean.orig_smi_canonical

    # The failing one falls back to the input molecule instead of propagating.
    assert len(salted.mols) == 1
    assert salted.mols[0] is salted.mol_orig_frm_inp_smi
    assert salted.orig_smi == "CCCCO.[Na+]"


def test_desalt_orig_smi_rejects_a_fragment_that_failed_sanitization(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Regression: the desalted molecule was added to the container without
    # checking that it survived sanitization, so a MyMol with rdkit_mol None
    # could reach the later steps and the writers.
    salted = MolContainer("CCCCO.[Na+]", "salted", 0, {})

    def unsanitizable_desalter(contnr: MolContainer) -> MyMol.MyMol:
        """Return a largest fragment that did not survive sanitization.

        Args:
            contnr: The container being desalted.

        Returns:
            A molecule whose rdkit_mol is None.
        """
        mol = MyMol.MyMol("CCCCO")
        mol.rdkit_mol = None
        return mol

    monkeypatch.setattr(DeSaltOrigSmiles, "desalter", unsanitizable_desalter)

    desalt_orig_smi([salted], 1, "serial", None)

    assert len(salted.mols) == 1
    kept = salted.mols[0]
    assert kept is salted.mol_orig_frm_inp_smi
    assert kept.rdkit_mol is not None
    # The fallback still has to stamp the provenance the skipped desalter would
    # have recorded.
    assert kept.genealogy == ["CCCCO.[Na+] (source)"]
