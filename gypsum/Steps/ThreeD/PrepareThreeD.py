"""Runs the 3D preparation process."""
import __future__

from gypsum.Steps.ThreeD.Convert2DTo3D import convert_2d_to_3d
from gypsum.Steps.ThreeD.GenerateAlternate3DNonaromaticRingConfs \
    import generate_alternate_3d_nonaromatic_ring_confs
from gypsum.Steps.ThreeD.Minimize3D import minimize_3d

def prepare_three_d(contnrs, params):
    """
    Runs the pipeline for generating the 3D models.
    """
    max_variants_per_compound = params["max_variants_per_compound"]
    thoroughness = params["thoroughness"]
    second_embed = params["second_embed"]

    num_processors = params["num_processors"]
    multithread_mode = params["multithread_mode"]
    Parallelizer_obj = params["Parallelizer"]

    if not params["2d_output_only"]:
        convert_2d_to_3d(contnrs, max_variants_per_compound, thoroughness, num_processors, multithread_mode, Parallelizer_obj)
    #self.print_current_smiles()

    if not params["skip_alternate_ring_conformations"]:
        generate_alternate_3d_nonaromatic_ring_confs(contnrs, thoroughness,
                                                     max_variants_per_compound,
                                                     num_processors,
                                                     second_embed, multithread_mode, Parallelizer_obj)
    #self.print_current_smiles()

    if not params["skip_optimize_geometry"]:
        minimize_3d(contnrs, thoroughness, max_variants_per_compound, num_processors, second_embed, multithread_mode, Parallelizer_obj)
    #self.print_current_smiles()
    