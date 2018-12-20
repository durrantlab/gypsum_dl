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
Runs the smile preparation process. Generates alternate ionization,
tautomeric, chiral forms, etc.
"""

import __future__

from gypsum import Utils
from gypsum.Steps.SMILES.DeSaltOrigSmiles import desalt_orig_smi
from gypsum.Steps.SMILES.AddHydrogens import add_hydrogens
from gypsum.Steps.SMILES.MakeTautomers import make_tauts
from gypsum.Steps.SMILES.EnumerateChiralMols import enumerate_chiral_molecules
from gypsum.Steps.SMILES.EnumerateDoubleBonds import enumerate_double_bonds

def prepare_smiles(contnrs, params):
    """Runs the appropriate steps for processing the SMILES strings.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param params: The user parameters.
    :type params: dict
    """

    # Unpack some of the parameter values.
    min_ph = params["min_ph"]
    max_ph = params["max_ph"]
    std_dev = params["ph_std_dev"]
    max_variants_per_compound = params["max_variants_per_compound"]
    thoroughness = params["thoroughness"]
    num_procs = params["num_processors"]
    job_manager = params["job_manager"]
    parallelizer_obj = params["Parallelizer"]

    debug = False

    # Desalt the molecules.
    Utils.log("Begin Desaltings")
    desalt_orig_smi(contnrs, num_procs, job_manager, parallelizer_obj)
    Utils.log("Done with Desalting")

    if debug: Utils.print_current_smiles(contnrs)

    # Add hydrogens for user-specified pH, if requested.
    if not params["skip_adding_hydrogen"]:
        Utils.log("Protonating Molecules")
        add_hydrogens(contnrs, min_ph, max_ph, std_dev, max_variants_per_compound,
                      thoroughness, num_procs, job_manager,
                      parallelizer_obj)
        Utils.log("Done with Protonating")
    else:
        Utils.log("Skipping Protonation")
        wrap_molecules(contnrs)

    if debug: Utils.print_current_smiles(contnrs)

    # Make alternate tautomeric forms, if requested.
    if not params["skip_making_tautomers"]:
        Utils.log("Tautomerizing Molecules")
        make_tauts(contnrs, max_variants_per_compound, thoroughness,
                   num_procs, job_manager, parallelizer_obj)
        Utils.log("Done with Tautomerization")
    else:
        Utils.log("Skipping Tautomerization")

    if debug: Utils.print_current_smiles(contnrs)

    # Make alternate chiral forms, if requested.
    if not params["skip_ennumerate_chiral_mol"]:
        Utils.log("Enumerating Chirality")
        enumerate_chiral_molecules(contnrs, max_variants_per_compound,
                                   thoroughness, num_procs,
                                   job_manager, parallelizer_obj)
        Utils.log("Done with Chirality Enumeration")
    else:
        Utils.log("Skipping Chirality Enumeration")

    if debug: Utils.print_current_smiles(contnrs)

    # Make alternate double-bond isomers, if requested.
    if not params["skip_ennumerate_double_bonds"]:
        Utils.log("Enumerating Double Bonds")
        enumerate_double_bonds(contnrs, max_variants_per_compound,
                               thoroughness, num_procs,
                               job_manager, parallelizer_obj)
        Utils.log("Done with Double Bond Enumeration")
    else:
        Utils.log("Skipping Double Bond Enumeration")

    if debug: Utils.print_current_smiles(contnrs)

def wrap_molecules(contnrs):
    """Each molecule container holds only one SMILES string
    (corresponding to the input structure). Dimorphite-DL can potentially
    produce multiple SMILES strings with different ionization states. There is
    no way to store muliple smiles in a molecule container.

    But those containers are designed to store multiple RDKit molecule
    objects. Gypsum stores Dimorphite-DL output molecules as RdKit molecule
    objects, in the container's mol list.

    But Gypsum gives the user the option of skipping the protonation step. In
    this case, the one SMILES needs to be converted to a RDKit mol object for
    subsequent steps to work. Let's do that here.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    """

    for mol_cont in contnrs:
        if len(mol_cont.mols) == 0:
            smi = mol_cont.orig_smi_canonical
            mol_cont.add_smiles(smi)
