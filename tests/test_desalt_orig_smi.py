"""Regression tests for SMILES desalting."""

from gypsum_dl.MolContainer import MolContainer
from gypsum_dl.steps.smiles.DeSaltOrigSmiles import desalt_orig_smi, desalter


def test_desalter_does_not_alias_original_genealogy() -> None:
    # Regression (bug 7): desalter shared the pristine molecule's genealogy list
    # (`new_mol.genealogy = contnr.mol_orig_frm_inp_smi.genealogy`), so the
    # later append in desalt_orig_smi mutated the container's original molecule.
    contnr = MolContainer("CCCCO.[Na+]", "salt", 0, {})
    orig = contnr.mol_orig_frm_inp_smi.genealogy
    before = list(orig)

    new_mol = desalter(contnr)
    assert new_mol.genealogy is not orig
    new_mol.genealogy.append("marker")
    assert contnr.mol_orig_frm_inp_smi.genealogy == before


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
