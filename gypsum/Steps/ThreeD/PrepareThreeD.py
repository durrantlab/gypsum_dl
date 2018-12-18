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
Runs the 3D preparation process.
"""

import __future__

from gypsum.Steps.ThreeD.Convert2DTo3D import convert_2d_to_3d
from gypsum.Steps.ThreeD.GenerateAlternate3DNonaromaticRingConfs \
    import generate_alternate_3d_nonaromatic_ring_confs
from gypsum.Steps.ThreeD.Minimize3D import minimize_3d

def prepare_3d(contnrs, params):
    """Runs the pipeline for generating the 3D small-molecule models.

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param params: The parameters.
    :type params: dict
    """

    # Unpack some of the user parameters.
    max_variants_per_compound = params["max_variants_per_compound"]
    thoroughness = params["thoroughness"]
    second_embed = params["second_embed"]
    num_procs = params["num_processors"]
    multithread_mode = params["multithread_mode"]
    parallelizer_obj = params["Parallelizer"]

    # Do the 2d to 3d conversionl, if requested.
    if not params["2d_output_only"]:
        # Make the 3D model.
        convert_2d_to_3d(contnrs, max_variants_per_compound, thoroughness,
                         num_procs, multithread_mode, parallelizer_obj)

        # Generate alternate non-aromatic ring conformations, if requested.
        if not params["skip_alternate_ring_conformations"]:
            generate_alternate_3d_nonaromatic_ring_confs(
                contnrs, max_variants_per_compound, thoroughness, num_procs,
                second_embed, multithread_mode, parallelizer_obj
            )

        # Minimize the molecules, if requested.
        if not params["skip_optimize_geometry"]:
            minimize_3d(contnrs, max_variants_per_compound, thoroughness, num_procs, second_embed, multithread_mode, parallelizer_obj)
