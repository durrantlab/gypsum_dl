"""Runs the 3D preparation process."""
from gypsum.Steps.ThreeD.Convert2DTo3D import convert_2d_to_3d
from gypsum.Steps.ThreeD.GenerateAlternate3DNonaromaticRingConfs \
    import generate_alternate_3d_nonaromatic_ring_confs
from gypsum.Steps.ThreeD.Minimize3D import minimize_3d

def prepare_three_d(self):
    """
    Runs the pipeline for generating the 3D models.
    """
    if not self.params["2d_output_only"]:
        convert_2d_to_3d(self)
    self.print_current_smiles()

    if not self.params["skip_alternate_ring_conformations"]:
        generate_alternate_3d_nonaromatic_ring_confs(self)
    self.print_current_smiles()

    if not self.params["skip_optimize_geometry"]:
        minimize_3d(self)
    self.print_current_smiles()
