"""
Runs the 3D preparation process.
"""

from typing import TYPE_CHECKING, Any

from gypsum_dl.steps.conf.convert import convert_2d_to_3d
from gypsum_dl.steps.conf.minimize import minimize_3d
from gypsum_dl.steps.conf.rings import generate_alternate_3d_nonaromatic_ring_confs

if TYPE_CHECKING:
    from gypsum_dl import MoleculeContainer


def prepare_3d(contnrs: list["MoleculeContainer"], params: dict[str, Any]) -> None:
    """Runs the pipeline for generating the 3D small-molecule models.

    Args:
        contnrs: A list of containers (container.MoleculeContainer).
        params: The parameters.
    """

    # Do the 2d to 3d conversion, if requested.
    if params["2d_output_only"]:
        return

    # Unpack some of the user parameters.
    max_variants_per_compound = params["max_variants_per_compound"]
    thoroughness = params["thoroughness"]
    num_procs = params["num_processors"]
    job_manager = params["job_manager"]
    parallelizer_obj = params["Parallelizer"]

    # Make the 3D model.
    convert_2d_to_3d(
        contnrs,
        max_variants_per_compound,
        thoroughness,
        num_procs,
        job_manager,
        parallelizer_obj,
    )

    second_embed = params["second_embed"]
    # Generate alternate non-aromatic ring conformations, if requested.
    if not params["skip_alternate_ring_conformations"]:
        generate_alternate_3d_nonaromatic_ring_confs(
            contnrs,
            max_variants_per_compound,
            thoroughness,
            num_procs,
            second_embed,
            job_manager,
            parallelizer_obj,
        )

    # Minimize the molecules, if requested.
    if not params["skip_optimize_geometry"]:
        minimize_3d(
            contnrs,
            max_variants_per_compound,
            thoroughness,
            num_procs,
            second_embed,
            job_manager,
            parallelizer_obj,
        )
