"""
This module generates alternate non-aromatic ring conformations on the fly,
since most modern docking programs (e.g., Vina) can't consider alternate ring
conformations.
"""

import __future__

import copy
import warnings
from collections.abc import Hashable

import gypsum_dl.parallelizer as Parallelizer
from gypsum_dl import chem_utils, utils
from gypsum_dl.MyMol import MyConformer, MyMol
from gypsum_dl.steps.conf.Minimize3D import parallel_minit

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except Exception:
    utils.exception("You need to install rdkit and its dependencies.")

try:
    import numpy
except Exception:
    utils.exception("You need to install numpy and its dependencies.")

try:
    from scipy.cluster.vq import kmeans2
except Exception:
    utils.exception("You need to install scipy and its dependencies.")


def generate_alternate_3d_nonaromatic_ring_confs(
    contnrs,
    max_variants_per_compound,
    thoroughness,
    num_procs,
    second_embed,
    job_manager,
    parallelizer_obj,
    minimize: bool = True,
) -> frozenset[int]:
    """Docking programs like Vina rotate chemical moieties around their
       rotatable bonds, so it's not necessary to generate a larger rotomer
       library for each molecule. The one exception to this rule is
       non-aromatic rings, which can assume multiple conformations (boat vs.
       chair, etc.). This function generates a few low-energy ring structures
       for each molecule with a non-aromatic ring(s).

    :param contnrs: A list of containers (MolContainer.MolContainer).
    :type contnrs: list
    :param max_variants_per_compound: To control the combinatorial explosion,
       only this number of variants (molecules) will be advanced to the next
       step.
    :type max_variants_per_compound: int
    :param thoroughness: How many molecules to generate per variant (molecule)
       retained, for evaluation. For example, perhaps you want to advance five
       molecules (max_variants_per_compound = 5). You could just generate five
       and advance them all. Or you could generate ten and advance the best
       five (so thoroughness = 2). Using thoroughness > 1 increases the
       computational expense, but it also increases the chances of finding good
       molecules.
    :type thoroughness: int
    :param num_procs: The number of processors to use.
    :type num_procs: int
    :param second_embed: Whether to try to generate 3D coordinates using an
        older algorithm if the better (default) algorithm fails. This can add
        run time, but sometimes converts certain molecules that would
        otherwise fail.
    :type second_embed: bool
    :param job_manager: The multiprocess mode.
    :type job_manager: string
    :param parallelizer_obj: The Parallelizer object.
    :type parallelizer_obj: Parallelizer.Parallelizer
    :param minimize: Whether to minimize the geometries of the ring conformers
        that are generated. Defaults to True.
    :type minimize: bool
    :return: The contnr_idx of every container with non-aromatic rings for
        which no ring conformer could be generated. Those containers keep the
        molecules they arrived with, which this step has not minimized.
    :rtype: frozenset[int]
    """

    # Let the user know you've started this step.
    utils.log(
        "Generating several conformers of molecules with non-aromatic "
        + "rings (boat vs. chair, etc.)..."
    )

    # A cap of zero means "do not enumerate variants," not "emit no models,"
    # so one ring conformer is still needed; a cap of zero would otherwise
    # reach kmeans2 with zero clusters.
    variant_cap = max(1, max_variants_per_compound)

    # Create parameters (inputs) to feed to the parallelizer.
    params = []
    ones_with_nonaro_rngs = set([])  # This is just to keep track of which
    # ones have non-aromatic rings.
    for contnr in contnrs:
        if contnr.num_nonaro_rngs > 0:
            # contnr_idx, not the list position: this set is later compared
            # against the keys of grouped, which come from mol.contnr_idx.
            ones_with_nonaro_rngs.add(contnr.contnr_idx)
            params.extend(
                (mol, variant_cap, thoroughness, second_embed, minimize)
                for mol in contnr.mols
            )
    params = tuple(params)

    # If there are no compounds with non-aromatic rings, no need to continue.
    if not ones_with_nonaro_rngs:
        return frozenset()  # There are no such ligands to process.

    # Run it through the parallelizer
    tmp = []
    if parallelizer_obj is None:
        # MultiThreading with one processor is the same dispatch serial mode
        # uses, and it draws a seed per job; calling run_one directly left
        # this path drawing from whatever generator state happened to be in
        # place, so --random_seed did not reach it.
        tmp = Parallelizer.MultiThreading(params, 1, parallel_get_ring_confs)
    else:
        tmp = parallelizer_obj.run(
            params, parallel_get_ring_confs, num_procs, job_manager
        )
    # Flatten the results. strip_none guards against tasks that returned None
    # (failed embedding); the worker also reports raised exceptions as None.
    results = Parallelizer.flatten_list(Parallelizer.strip_none(tmp))

    # Group by mol. You can't use existing functions because they would
    # require you to recalculate already calculated energies.
    grouped = {}  # Index will be container index. Value is list of
    # (energy, mol) pairs.
    for mol in results:
        # Save the energy as a prop while you're here. The raw value, sentinel
        # included, is what the sort below needs; only the published property
        # drops it.
        energy = mol.conformers[0].energy
        mol.mol_props["Energy"] = utils.energy_for_output(energy)

        # Add the mol with it's energy to the appropriate entry in grouped.
        # Make that entry if needed.
        contnr_idx = mol.contnr_idx
        if contnr_idx not in grouped:
            grouped[contnr_idx] = []
        grouped[contnr_idx].append((energy, mol))

    # Now, for each container, keep only the best ones.
    contnr_by_idx = utils.contnrs_by_idx(contnrs)
    for contnr_idx, lst_enrgy_mol_pairs in grouped.items():
        contnr = contnr_by_idx[contnr_idx]
        prior_mols = contnr.mols
        contnr.mols = []  # Note that only affects ones that
        # had non-aromatic rings.
        for mol in _select_ring_confs_by_form(
            lst_enrgy_mol_pairs, prior_mols, variant_cap
        ):
            contnr.add_mol(mol)

    # Any container that had non-aromatic rings but produced no results (all
    # ring-conformer generation failed) is absent from grouped. Its original
    # mols are untouched; flag them so the failure is recorded.
    failed_contnr_idxs = frozenset(ones_with_nonaro_rngs - set(grouped.keys()))
    for contnr_idx in failed_contnr_idxs:
        for mol in contnr_by_idx[contnr_idx].mols:
            mol.genealogy.append(
                "(WARNING: Could not generate alternate conformations "
                + "of nonaromatic ring)"
            )

    return failed_contnr_idxs


def _form_key(mol: MyMol) -> Hashable:
    """Identify which chemical form (protonation state, tautomer, isomer) a
    ring conformer belongs to.

    Args:
        mol: A ring conformer or an input variant.

    Returns:
        The canonical SMILES, or the molecule itself when that could not be
            computed, so such a molecule is a form of its own (matching the
            identity fallback in MyMol.__eq__ and __hash__).
    """
    smi = mol.smiles()
    return smi if isinstance(smi, str) else mol


def _select_ring_confs_by_form(
    pairs: list[tuple[float, MyMol]], prior_mols: list[MyMol], cap: int
) -> list[MyMol]:
    """Pick up to cap ring conformers for one container, spreading the picks
    across chemical forms.

    UFF energies are comparable only between conformers of the same form.
    Protonation states and tautomers differ in atoms or bonding, so ranking
    them against each other let an energy offset with no physical meaning
    decide which forms survived, overriding the choices the SMILES steps had
    made. Here each form's conformers are ranked by energy among themselves,
    and the forms take turns: every form's best conformer first, then every
    form's second best, and so on until the cap is reached.

    Args:
        pairs: (energy, molecule) for every ring conformer of the container.
        prior_mols: The container's molecules as they entered this step. Forms
            take their turns in this order, which is the order the earlier
            steps left them in. Forms not found there follow, by SMILES, so
            the result never depends on hash or parallelizer ordering.
        cap: The largest number of conformers to keep.

    Returns:
        The selected molecules, in the order they were picked.
    """
    by_form: dict[Hashable, list[tuple[float, MyMol]]] = {}
    for pair in pairs:
        by_form.setdefault(_form_key(pair[1]), []).append(pair)

    rank: dict[Hashable, int] = {}
    for i, mol in enumerate(prior_mols):
        rank.setdefault(_form_key(mol), i)
    forms = sorted(
        by_form,
        key=lambda form: (
            rank.get(form, len(prior_mols)),
            by_form[form][0][1].smiles() or "",
        ),
    )

    # Conformers of one form share a SMILES, so energy alone orders them;
    # the stable sort leaves exact ties in the order the worker produced.
    queues = [sorted(by_form[form], key=lambda pair: pair[0]) for form in forms]

    selected: list[MyMol] = []
    depth = 0
    while len(selected) < cap and any(depth < len(queue) for queue in queues):
        for queue in queues:
            if depth < len(queue) and len(selected) < cap:
                selected.append(queue[depth][1])
        depth += 1
    return selected


def parallel_get_ring_confs(
    mol, max_variants_per_compound, thoroughness, second_embed, minimize: bool = True
):
    """Gets alternate ring conformations. Meant to run with the parallelizer class.

    :param mol: The molecule to process (with non-aromatic ring(s)).
    :type mol: MyMol.MyMol
    :param max_variants_per_compound: To control the combinatorial explosion,
       only this number of variants (molecules) will be advanced to the next
       step.
    :type max_variants_per_compound: int
    :param thoroughness: How many molecules to generate per variant (molecule)
       retained, for evaluation. For example, perhaps you want to advance five
       molecules (max_variants_per_compound = 5). You could just generate five
       and advance them all. Or you could generate ten and advance the best
       five (so thoroughness = 2). Using thoroughness > 1 increases the
       computational expense, but it also increases the chances of finding good
       molecules.
    :type thoroughness: int
    :param second_embed: Whether to try to generate 3D coordinates using an
        older algorithm if the better (default) algorithm fails. This can add
        run time, but sometimes converts certain molecules that would
        otherwise fail.
    :type second_embed: bool
    :param minimize: Whether to minimize the geometries of the generated
        conformers. Defaults to True.
    :type minimize: bool
    :return: A list of MyMol.MyMol objects, with alternate ring conformations.
    :rtype: list
    """

    # Serial and in-process runs pass the container's own molecule, while
    # multiprocessing passes a pickled copy. Work on a copy in every mode, so
    # a molecule whose conformer search fails part way is carried over in the
    # same state regardless of the job manager.
    mol = copy.deepcopy(mol)

    # Make it easier to access the container index.
    contnr_idx = mol.contnr_idx

    # All the molecules in this container must have nonatomatic rings (because
    # they are all variants of the same source molecule). So just make a new
    # mols list.

    # Get the ring atom indecies
    rings = mol.get_idxs_of_nonaro_rng_atms()

    # It's possible a variant (e.g., a tautomer) is aromatic even if the
    # original molecule was not. In that case, there are no non-aromatic
    # rings to generate conformers for.
    if len(rings) == 0:
        if not minimize:
            return [mol]
        # minimize_3d leaves every molecule in a ring-bearing container to
        # this step, so a variant with no ring to vary still has to be
        # minimized here or it ships with its raw embedded geometry.
        minimized = parallel_minit(
            mol, max_variants_per_compound, thoroughness, second_embed
        )
        return None if minimized is None else [minimized]

    # Convert that into the bond indecies.

    # A list of lists, where each inner list has the indexes of the bonds that
    # comprise a ring.
    rings_by_bond_indexes = []
    for ring_atom_indecies in rings:
        bond_indexes = []
        for ring_atm_idx in ring_atom_indecies:
            a = mol.rdkit_mol.GetAtomWithIdx(ring_atm_idx)
            bonds = a.GetBonds()
            for bond in bonds:
                atom_indecies = [bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]
                atom_indecies.remove(ring_atm_idx)
                other_atm_idx = atom_indecies[0]
                if other_atm_idx in ring_atom_indecies:
                    bond_indexes.append(bond.GetIdx())
        bond_indexes = sorted(set(bond_indexes))
        rings_by_bond_indexes.append(bond_indexes)

    # Generate a bunch of conformations, ordered from best energy to worst.
    # Note that this is cached. Minimizing too, unless the caller asked for the
    # optimization step to be skipped.
    mol.add_conformers(thoroughness * max_variants_per_compound, 0.1, minimize)

    if len(mol.conformers) > 0:
        # Sometimes there are no conformers if it's an impossible structure.
        # Like
        # [H]c1nc(N2C(=O)[C@@]3(C([H])([H])[H])[C@@]4([H])O[C@@]([H])(C([H])([H])C4([H])[H])[C@]3(C([H])([H])[H])C2=O)sc1[H]
        # So don't save this one anyway.

        # Get the scores (lowest energy) of these minimized conformers.
        mol.load_conformers_into_rdkit_mol()

        # Extract just the rings.
        ring_mols = [
            Chem.PathToSubmol(mol.rdkit_mol, bi) for bi in rings_by_bond_indexes
        ]

        # Align get the rmsds relative to the first conformation, for each
        # ring separately.
        list_of_rmslists = [[]] * len(ring_mols)
        for k in range(len(ring_mols)):
            list_of_rmslists[k] = []
            AllChem.AlignMolConformers(ring_mols[k], RMSlist=list_of_rmslists[k])

        # Get points for each conformer (rmsd_ring1, rmsd_ring2, rmsd_ring3)
        pts = numpy.array(list_of_rmslists).T

        # With no ring submols to align, list_of_rmslists is empty and the
        # transpose has no second dimension, so there is nothing to cluster on.
        # A single conformer is not this case: every RMS list is then empty and
        # pts has shape (0, num_rings), which the vstack below turns into the
        # one row describing the reference conformer.
        if pts.ndim != 2 or pts.shape[1] == 0:
            return [mol]

        pts = numpy.vstack((numpy.array([[0.0] * pts.shape[1]]), pts))
        # Cluster those points, get lowest-energy member of each.
        if len(pts) < max_variants_per_compound:
            num_clusters = len(pts)
        else:
            num_clusters = max_variants_per_compound

        # When kmeans2 runs on insufficient clusters, it can sometimes throw an
        # error about empty clusters. This is not necessary to throw for the
        # user and so we have supressed it here.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            groups = kmeans2(pts, num_clusters, minit="points")[1]

        # These are geometrically diverse conformations of this one form. The
        # other forms in the container (enantiomers, tautomers, etc.) are
        # clustered separately, and the final selection takes conformers from
        # each form in turn.

        # Key is group id from kmeans (int). Values are the MyMol.MyConformers
        # objects.
        best_conf_per_group = {}

        conformers = mol.rdkit_mol.GetConformers()
        for k, grp in enumerate(groups):
            if grp not in list(best_conf_per_group.keys()):
                best_conf_per_group[grp] = mol.conformers[k]
        # best_confs has the MyMol.MyConformers objects.
        best_confs = best_conf_per_group.values()

        # Convert rdkit mols to MyMol.MyMol and save those MyMol.MyMol objects
        # for returning.
        results = []
        for conf in best_confs:
            new_mol = copy.deepcopy(mol)
            c = MyConformer(new_mol, conf.conformer(), second_embed)
            new_mol.conformers = [c]
            energy = c.energy

            new_mol.genealogy = mol.genealogy[:]
            new_mol.genealogy.append(
                new_mol.smiles(True)
                + " (nonaromatic ring conformer: "
                + utils.describe_energy(energy)
                + ")"
            )

            results.append(new_mol)  # i is mol index

        return results

    # If you get here, something went wrong.
    return None
