"""
This module is for testing Gypsum-DL. Not quite unit tests, but good enough
for now.
"""

import glob
import os
import shutil

from gypsum_dl.start import prepare_molecules


def test_samples(test_dir):
    path_smiles = os.path.join(test_dir, "files/sample/sample_molecules.smi")
    output_folder = os.path.join(test_dir, "tmp/sample")

    # Delete test output directory if it exists.
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    # Make the directory
    os.mkdir(output_folder)

    # Make the Gypsum-DL parameters.
    params = {
        "source": path_smiles,
        "separate_output_files": True,
        "job_manager": "serial",  # multiprocessing
        "output_folder": output_folder,
        "add_pdb_output": False,
        "max_variants_per_compound": 8,
        "thoroughness": 1,
        "min_ph": 4,
        "max_ph": 10,
        "pka_precision": 1,
        "use_durrant_lab_filters": True,
    }

    # Prepare the molecules.
    prepare_molecules(params)

    # Get the output sdf files.
    sdf_files = glob.glob(f"{output_folder}/*")

    assert len(sdf_files) == 15, f"Expected 15 output files, got {len(sdf_files)}."

    # Get all the smiles from the files.
    all_smiles = set([])
    for sdf_file in sdf_files:
        lines = open(sdf_file).readlines()
        for i, line in enumerate(lines):
            if "<SMILES>" in line:
                all_smiles.add(lines[i + 1].strip())

    # List what the smiles should be.
    target_smiles = set([])

    # salt_and_ionization should produce two models (ionized and
    # deionized).
    target_smiles |= {"[O-]c1ccccc1", "Oc1ccccc1"}

    # tautomer_and_cis_trans should produce three models (two tautomers, one
    # of them with alternate cis/trans).
    target_smiles |= {r"C/C=C\O", "C/C=C/O", "CCC=O"}

    # two_chiral_one_unspecified_and_tautomer should produce four models.
    target_smiles |= {
        "CC(C)C(=O)[C@@](F)(Cl)C[C@@](C)(F)Cl",
        "CC(C)=C(O)[C@@](F)(Cl)C[C@@](C)(F)Cl",
        "CC(C)C(=O)[C@](F)(Cl)C[C@@](C)(F)Cl",
        "CC(C)=C(O)[C@](F)(Cl)C[C@@](C)(F)Cl",
    }

    # two_double_bonds_one_chiral_center should produce eight models.
    target_smiles |= {
        r"CC/C(C[C@@](C)(Cl)I)=C(I)\C(F)=C(/C)Cl",
        "CC/C(C[C@](C)(Cl)I)=C(I)/C(F)=C(/C)Cl",
        r"CC/C(C[C@](C)(Cl)I)=C(I)/C(F)=C(\C)Cl",
        r"CC/C(C[C@](C)(Cl)I)=C(I)\C(F)=C(\C)Cl",
        r"CC/C(C[C@@](C)(Cl)I)=C(I)/C(F)=C(\C)Cl",
        r"CC/C(C[C@@](C)(Cl)I)=C(I)\C(F)=C(\C)Cl",
        "CC/C(C[C@@](C)(Cl)I)=C(I)/C(F)=C(/C)Cl",
        r"CC/C(C[C@](C)(Cl)I)=C(I)\C(F)=C(/C)Cl",
    }

    # two_double_bonds_one_unspecified should produce two models.
    target_smiles |= {
        r"CC/C(C)=C(\Cl)C/C(I)=C(\C)F",
        r"CC/C(C)=C(/Cl)C/C(I)=C(\C)F",
    }

    # non_aromatic_ring should produce one model. It will list it several
    # times, because different ring conformations of the same model.
    target_smiles |= {"CC(C)(C)[C@H]1CC[C@@H](C(C)(C)C)CC1"}

    # There should be no =[N-] if Durrant lab filters are turned on. Note:
    # Removed "CC(=N)O" from below list because durrant lab filters now remove
    # iminols.
    target_smiles |= {"CC([NH-])=O", "CC(N)=O"}

    # There should be no [N-]C=[N+] (CC(=O)[N-]C=[N+](C)C).
    target_smiles |= {
        r"CC(=O)/N=C\[NH+](C)C",
        "CC(=O)/N=C/[NH+](C)C",
        "CC(=O)NC=[N+](C)C",
    }

    # There should be no [nH+]c[n-] (c1c[nH+]c[n-]1)
    target_smiles |= {"c1c[n-]cn1", "c1c[nH+]c[nH]1", "c1c[nH]cn1"}

    # There should be no [#7+]~[#7+] (c1cc[nH+][nH+]c1)
    target_smiles |= {"c1ccnnc1", "c1cc[nH+]nc1"}

    # There should be no [#7-]~[#7-] (CC(=O)[N-][N-]C(C)=O). Note that some
    # are commented out because Python2 and Python3 given different SMILES
    # strings that are all valid. See below to see how things are
    # consolodated. (Really this was probably a bad example to pick because
    # there are so many forms...)
    target_smiles |= {"CC(=O)NNC(C)=O"}
    # r"CC(=O)N/N=C(\C)O",
    # r"CC(=O)[N-]/N=C(/C)O",
    # r"C/C(O)=N/N=C(\C)O",
    # r"C/C(O)=N\N=C(/C)O",  # No longer allowing internal iminols
    # r"CC(=O)[N-]/N=C(\C)O",
    # "CC(=O)[N-]NC(C)=O",
    # "CC(=O)N/N=C(/C)O"

    # There should be no [!#7]~[#7+]~[#7-]~[!#7] (c1c[n-][nH+]c1)
    target_smiles |= {"c1cn[n-]c1", "c1cn[nH]c1", "c1c[nH][nH+]c1"}

    # Azides can have adjacent +/- nitrogens.
    target_smiles |= {"CN=[N+]=[N-]", "CN=[N+]=N"}

    # Python3 gives some smiles that are different than thsoe obtain with
    # Python2. But they are just different representations of the same thing.
    # Let's make the switch to the Python2 form for this test.
    all_smiles = {"CN=[N+]=N" if s == "[H]N=[N+]=NC" else s for s in all_smiles}

    # Note: Commented out below because durrant lab filters now remove
    # iminols.
    # all_smiles = set(
    #     ["CC(=N)O" if s in [r"[H]/N=C(\C)O", "[H]/N=C(/C)O"] else s for s in all_smiles]
    # )

    all_smiles = {
        # Different one that turns up sometimes
        r"C/C(O)=N\N=C(/C)O" if s == r"C/C(O)=N/N=C(/C)O" else s
        for s in all_smiles
    }
    all_smiles = {
        (
            r"CC(=O)NNC(C)=O"
            if s
            in [
                r"CC(=O)[N-]/N=C(\C)O",
                r"C/C(O)=N/N=C(\C)O",
                r"CC(=O)N/N=C(\C)O",
                r"CC(=O)[N-]/N=C(/C)O",
                r"CC(=O)[N-]NC(C)=O",
                r"CC(=O)N/N=C(/C)O",
            ]
            # Different one that turns up sometimes
            else s
        )
        for s in all_smiles
    }

    assert len(all_smiles) == len(target_smiles)

    assert len(all_smiles ^ target_smiles) == 0, (
        f"Differences in smiles: {list(all_smiles ^ target_smiles)}"
    )
