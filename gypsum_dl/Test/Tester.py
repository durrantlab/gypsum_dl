# Copyright 2018 Jacob D. Durrant
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
This module is for testing Gypsum-DL. Not quite unit tests, but good enough
for now.
"""

import os
import shutil
import glob
from gypsum_dl import Utils
from gypsum_dl.Start import prepare_molecules

def run_test():
    script_dir = os.path.dirname(os.path.realpath(__file__))
    output_folder = script_dir + os.sep + "gypsum_dl_test_output" + os.sep

    # Delete test output directory if it exists.
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)

    # Make the directory
    os.mkdir(output_folder)

    # Make the Gypsum-DL parameters.
    params = {
        "source": script_dir + os.sep + "sample_molecules.smi",
        "separate_output_files": True,
        "job_manager": "serial",  # multiprocessing
        "output_folder": output_folder,
        "add_pdb_output": False,
        "max_variants_per_compound": 8,
        "thoroughness": 1,
        "min_ph": 4,
        "max_ph": 10,
        "pka_precision": 1,
        "use_durrant_lab_filters": True
    }

    # Prepare the molecules.
    prepare_molecules(params)
    Utils.log("")
    Utils.log("TEST RESULTS")
    Utils.log("============")

    # Get the output sdf files.
    sdf_files = glob.glob(output_folder + "*")

    # There should be seven sdf files.
    msg = "Expected 7 output files, got " + str(len(sdf_files)) + "."
    if len(sdf_files) != 7:
        Utils.exception("FAILED. " + msg)
    else:
        Utils.log("PASSED. " + msg)

    # Get all the smiles from the files.
    all_smiles = set([])
    for sdf_file in sdf_files:
        lines = open(sdf_file).readlines()
        for i, line in enumerate(lines):
            if "<SMILES>" in line:
                all_smiles.add(lines[i+1].strip())

    # List what the smiles should be.
    target_smiles = set([])

    # salt_and_ionization should produce two models (ionized and
    # deionized).
    target_smiles |= set(["[O-]c1ccccc1", "Oc1ccccc1"])

    # tautomer_and_cis_trans should produce three models (two tautomers, one
    # of them with alternate cis/trans).
    target_smiles |= set(["C/C=C\O", "C/C=C/O", "CCC=O"])

    # two_chiral_one_unspecified_and_tautomer should produce four models.
    target_smiles |= set([
        "CC(C)C(=O)[C@@](F)(Cl)C[C@@](C)(F)Cl",
        "CC(C)=C(O)[C@@](F)(Cl)C[C@@](C)(F)Cl",
        "CC(C)C(=O)[C@](F)(Cl)C[C@@](C)(F)Cl",
        "CC(C)=C(O)[C@](F)(Cl)C[C@@](C)(F)Cl"
    ])

    # two_double_bonds_one_chiral_center should produce eight models.
    target_smiles |= set([
        "CC/C(C[C@@](C)(Cl)I)=C(I)\C(F)=C(/C)Cl",
        "CC/C(C[C@](C)(Cl)I)=C(I)/C(F)=C(/C)Cl",
        "CC/C(C[C@](C)(Cl)I)=C(I)/C(F)=C(\C)Cl",
        "CC/C(C[C@](C)(Cl)I)=C(I)\C(F)=C(\C)Cl",
        "CC/C(C[C@@](C)(Cl)I)=C(I)/C(F)=C(\C)Cl",
        "CC/C(C[C@@](C)(Cl)I)=C(I)\C(F)=C(\C)Cl",
        "CC/C(C[C@@](C)(Cl)I)=C(I)/C(F)=C(/C)Cl",
        "CC/C(C[C@](C)(Cl)I)=C(I)\C(F)=C(/C)Cl"
    ])

    # two_double_bonds_one_unspecified should produce two models.
    target_smiles |= set([
        "CC/C(C)=C(\Cl)C/C(I)=C(\C)F", "CC/C(C)=C(/Cl)C/C(I)=C(\C)F"
    ])

    # non_aromatic_ring should produce one model. It will list it several
    # times, because different ring conformations of the same model.
    target_smiles |= set([
        "CC(C)(C)[C@H]1CC[C@@H](C(C)(C)C)CC1"
    ])

    # msg = "Expected " + str(len(target_smiles)) + " total SMILES, got " + \
    #     str(len(all_smiles)) + "."
    # if len(all_smiles) != len(target_smiles):
    #     Utils.exception("FAILED. " + msg)
    # else:
    #     Utils.log("PASSED. " + msg)

    if len(all_smiles ^ target_smiles) > 0:
        Utils.exception("FAILED. " + "Got some SMILES I didn't expect: " + \
            " ".join(list(all_smiles ^ target_smiles)))
    else:
        Utils.log("PASSED. Gypsum-DL output the very SMILES strings I was expecting.")

    Utils.log("")

    # Delete test output directory if it exists.
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
