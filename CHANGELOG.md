# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## WIP: [2.0.0] - 2026-09-30

This release focuses on correctness, consistency, and reproducibility. Some of
these improvements affect which variants Gypsum-DL generates, so output may
differ from 2.0.0, including for runs that use the same random seed.

### Added

-   Added the `--random_seed` parameter. Seeds are drawn per job in the
    dispatching process, so a seeded run produces the same output under the
    serial, multiprocessing, and MPI job managers.
-   Added a `python -m gypsum_dl` entry point. MPI runs now launch with
    `mpirun -n N python -m mpi4py -m gypsum_dl ...`.
-   In MPI mode, the root process now writes `gypsum_dl_failed__inputN.smi` for
    jobs that raised, matching the other job managers.
-   Added validation for `job_manager`, `pH`, `thoroughness`,
    `max_variants_per_compound`, and the input source. `thoroughness` must now
    be between 1 and 1000.

### Changed

-   Seeded output differs from earlier releases because per-job seeding and
    cached ranking energies change the sequence of random draws.
-   Ring conformers are now selected within each chemical form (protonation
    state or tautomer) rather than by pooled UFF energy, since UFF energies are
    not comparable across forms. Every form kept by the SMILES steps now gets
    at least one output slot.
-   For molecules without non-aromatic rings, every embedded conformer is now
    minimized before the lowest-energy one is kept. Previously only the top
    `max_variants_per_compound` conformers by unminimized UFF energy were
    minimized, but that energy does not predict the minimized ranking, so
    `thoroughness` had little effect when `max_variants_per_compound` was
    small. The final minimization step now takes roughly 30 to 45 percent
    longer.
-   Ionization states are now chosen by Dimorphite-DL's probability ranking.
    Gypsum-DL previously requested `thoroughness × max_variants_per_compound`
    states and kept the ones with the lowest unminimized UFF energy, which
    cannot compare states that differ in atom count and charge. It now
    requests `max_variants_per_compound` states, so `thoroughness` no longer
    affects ionization. This requires Dimorphite-DL 2.1.0 or later.
-   Variants are now ranked by UFF energy only against variants with the same
    molecular formula and net charge, and each such group (in practice, each
    protonation state) gets a slot before any group gets a second. The
    tautomer, enantiomer, cis-trans, and first 3D culls previously ranked all
    of a compound's variants together, so a protonation state could be dropped
    because its energy, which is not comparable across states, was higher.
-   The Durrant-lab metal filter now matches atomic numbers on the variant
    structure rather than SMILES substrings. This keeps krypton from matching
    `[K` and lets the filter recognize isotope-labeled metals.
-   Durrant-lab SMARTS queries are now compiled once per process rather than
    being sent with each container.
-   The tautomer enumeration budget now scales with `thoroughness`.
-   Chiral-center detection and cis/trans enumeration now use RDKit's legacy
    stereo perception, so results are consistent across RDKit versions.
-   The PDB "Final SMILES" remark now matches the SMILES written to the SDF and
    HTML outputs.
-   Fields that Gypsum-DL computes now take precedence over tags in input SDF
    files.
-   HTML output is now written as UTF-8.
-   mpi4py is now optional. It is only needed for MPI mode.
-   `max_variants_per_compound = 0` is now supported and documented: Gypsum-DL
    desalts, skips SMILES enumeration, and outputs one 3D model per compound.
-   Updated the help text, error messages, `README.md`, and Pitt CRC example
    scripts to reflect the new MPI launch command. The documentation no longer
    describes multiprocessing and MPI runs as nondeterministic.
-   Warnings about steps that produced no variants now name what the step was
    generating (tautomers, enantiomers, etc.).

### Removed

-   Removed the unused `cache_prerun` parameter. Passing it now raises an
    unrecognized-parameter error.
-   Removed `end_time` and `run_time` from the saved parameters, since the
    parameter record is written before the run finishes. Both are still
    reported in the log.
-   Removed an unused tautomer filter, `durrant_lab_contains_bad_substr`,
    `prohibited_smi_substrs_for_substr`, and several unused `MolObjectHandling`
    helpers.

### Fixed

-   Parallel execution and error handling:
    -   Multiprocessing runs now continue when a job raises, and Gypsum-DL
        detects workers that exit unexpectedly.
    -   Serial, multiprocessing, and MPI jobs now share one failure handler,
        so an error in one job is isolated from the rest, and every failed job
        is reported.
    -   Improved MPI per-rank failure files, result reindexing, and handling
        of the parameters SDF.
    -   Enumeration, 3D, and filter workers now leave their inputs unchanged,
        so serial and multiprocessing runs give the same results.
-   Variant enumeration and selection:
    -   Corrected list indexing during low-energy variant selection.
    -   Containers are now indexed by `contnr_idx` rather than list position.
        Input compounds that share a SMILES string are now processed
        independently, and deduplication stays within each compound.
    -   Steps that generate no variants now consistently keep the existing
        structures, and tautomer-step failures are now recorded in the
        genealogy.
    -   Removed a redundant metal pre-filter in favor of the full Durrant-lab
        filter step, which already removes these compounds.
    -   Improved chiral-center counting, handling of truncated chirality
        assignments, and the chirality reference used by the tautomer filter.
    -   Improved double-bond enumeration for molecules with only one possible
        combination, deduplication in molecules without explicit hydrogens,
        the terminal-alkene stereo check, the enumeration budget, and
        reproducibility.
    -   Improved desalting tie-breaks, fragment pairing, and genealogy
        records. Desalting is now deterministic.
    -   Improved ionization fallbacks and provenance tracking, and preserved
        names in `add_smiles`. The genealogy now records the source SMILES,
        and each variant's `UniqueID` is now unique.
    -   SMILES caches are now refreshed when a molecule changes, including when
        hydrogen handling renumbers atoms.
-   3D generation and minimization:
    -   `second_embed` and `skip_optimize_geometry` are now applied as
        documented.
    -   Ring-containing molecules are now minimized when the ring-conformer
        step is skipped. Ringless variants, and molecules for which ring
        conformers could not be generated, are now minimized as well.
    -   Improved conformer ordering after minimization, conformer IDs, the
        ring-conformer shape check, and handling of missing conformers in
        RMSD calculation, conformer generation, and `minimize_3d`.
-   Input and output:
    -   Improved creation of nested `output_folder` paths, 2D SDF output,
        per-input SDF files, handling of atomless SDF records and `None` SDF
        properties, and handling of duplicate molecule names.
    -   Improved failure-file names, file encodings, `source_dir` handling,
        and file-handle cleanup. Output writers now skip molecules they
        cannot write.
-   Command line and logging:
    -   All command-line flags are now passed through correctly. Also
        improved boolean parameter checks, the default number of processors,
        and CLI documentation.
    -   Improved log environment variable parsing, line wrapping, newline
        handling, and debug-output settings. Sanitization failures are now
        reported without a SMILES string.
-   Error handling for problematic molecules:
    -   Improved handling of edge cases in `standardize_smiles`, `MyMol` SMILES
        generation, non-string canonical SMILES, and `flatten_list` results
        from failed workers.

## [1.3.0] - 2025-11-17

-   Bumping version to 1.3.0 to match the version on PyPI.

## [1.2.3] - 2025-11-05

-   Improved error handling for problematic input SMILES.

## [1.2.2] - 2025-11-03

### Changed

-   Refactored codebase to use pixi as our development environment and make this package pip installable.

### Fixed

-   Fixed an `IndexError` that rarely occurred during non-aromatic ring conformation generation.

## [1.2.1]

-   Fixed a bug in generating multiple non-aromatic ring conformations. This
  functionality worked correctly when Gypsum-DL was originally published (e.g.,
  with rdkit 2020.03.1). But due to changes made to more recent versions of
  rdkit, Gypsum-DL produced only one ring conformation, even if multiple were
  reasonably possible. It now produces multiple ring conformations even on
  recent versions of rdkit (e.g., `2023.03.1`). We recommend using the lastest
  version of rdkit.
-   Gypsum-DL now uses `AllChem.ETKDGv3` if it's available.
-   Modernized codebase some. Python2 no longer officially supported.
-   Updated copyright year to 2023.

## [1.2.0]

-   Added to the Durrant-lab filters to compensate for an amide-related bug in
  MolVS, one of Gypsum-DL's dependencies. MolVS sometimes tautomerizes
  `NC(=O)C[*]` to `N\C(O)=C\[*]`, so the Durrant-lab filters now remove any
  tautomers with substructures that match the SMARTS string
  `[$(N)]C(=C)[$([OH]),$([O-])]`.
-   Previously, the Durrant-lab filters only removed terminal iminols, which are
  improbable tautomers of terminal amides. According to
  [DataWarrior](https://openmolecules.org/datawarrior/), internal iminols (e.g.,
  `C\N=C(\C)O`) are also improbable, so these are now removed as well.

## [1.1.9]

-   Improved error handling when loading SDF files that are poorly formatted
  (e.g., that do not specify charged nitrogen atoms). Gypsum-DL depends on
  RDKit for SDF loading, and RDKit apparently cannot handle these errors. If
  you find that Gypsum-DL skips many of your compounds with a `Warning: Could
  not convert some SDF-formatted files to SMILES...` error, consider using an
  SMI (SMILES) file instead.

## [1.1.8]

Updated `README.md` to help some users who were having trouble installing
RDKit.

## [1.1.7]

-   Updated the `README.md` file, specifically the `Important Caveats` section.
-   Modest speed improvements when enumerating compounds with many chiral
  centers. (No need to enumerate far more compounds than will ultimately be
  used, given the values of the `thoroughness` and `max_variants_per_compound`
  user parameters.) This update should also allow Gypsum-DL to more
  efficiently use available memory.
-   Similar speed and memory improvements when enumerating compounds with many
  double bonds that have unspecified stereochemistries.

## [1.1.6]

-   Corrected minor bug that caused Durrant-lab filters to inappropriately
  retain some compounds when running in multiprocessing mode.
-   Fixed testing scripts, now that Durrant-lab filters remove iminols.

## [1.1.5]

-   Updated Dimorphite-DL to 1.2.4. Now better handles compounds with
  polyphosphate chains (e.g., ATP).
-   Minor updates to the Durrant-lab filters:
-   When running Gypsum-DL without the `--use_durrant_lab_filters` parameter,
    Gypsum-DL now displays a warning. We strongly recommend using these
    filters, but we choose not to turn them on by default in order to maintain
    backwards compatibility.
-   Added filter to compensate for a phosphate-related bug in MolVS, one of
    Gypsum-DL's dependencies. MolVS sometimes tautomerizes `[O]P(O)([O])=O` to
    `[O][PH](=O)([O])=O`, so the Durrant-lab filters now remove any tautomers
    with substructures that match the SMARTS string `O=[PH](=O)([#8])([#8])`.
-   Added filters to compensate for frequently seen, unusual MolVS
    tautomerization of adenine. The Durrant-lab filters now remove tautomers
    with substructures that match `[#7]=C1[#7]=C[#7]C=C1` and
    `N=c1cc[#7]c[#7]1`.
-   Added filter to remove terminal iminols. While amide-iminol
    tautomerization is valid, amides are far more common, and accounting for
    this tautomerization produces many improbable iminol compounds. The
    Durrant-lab filters now remove compounds with substructures that match
    `[$([NX2H1]),$([NX3H2])]=C[$([OH]),$([O-])]`.
-   Added filter to remove molecules containing `[Bi]`.
-   Gypsum-DL now outputs molecules with total charges between -4e and +4e.
  Before, the cutoff was -2e to 2e. We expanded the range to permit ATP and
  other similar molecules.

## [1.1.4]

-   Updated Dimorphite-DL to 1.2.3.
-   Added `sys.stdout.flush()` commands to ParallelMPI.run (see
  `gypsum_dl/gypsum_dl/Parallelizer.py`) to ensure that print statements
  properly output in large MPI runs.

## [1.1.3]

-   Gypsum-DL used to crash when provided with certain mal-formed SMILES
  strings. It now just skips those SMILES and warns the user that they are
  poorly formed. See Start.py:303 and MyMol.py:747.
-   Durrant-lab filters now remove molecules containing metal and boron atoms.
-   Some Durrant-lab filters are now applied immediately after desalting. We
  discovered that certain substructures cause Gypsum-DL to delay during the
  add-hydrogens step, specifically when Gypsum-DL generates the 3D structures
  required to rank conformers. Removing these compounds before adding
  hydrogens avoids the problem.
-   Improved code formatting.
-   Made minor spelling corrections to the output.

## [1.1.2]

-   Bug fix: thoroughness parameter is now correctly recognized as an int when
  specified from the command line.

## [1.1.1]

-   Updated Dimorphite-DL to version 1.2.2.
-   Corrected spelling in user-parameter names. Parameters that previously used
  "ennumerate" now use "enumerate".
-   Updated MolVS-generated tautomer filters. Previous versions of Gypsum-DL
  rejected tautomers that changed the number of _specified_ chiral centers. By
  default, Gypsum now rejects tautomers that change the total number of chiral
  centers, _both specified and unspecified_. To override the new default
  behavior (i.e., to allow tautomers that change the total number of chiral
  centers), use `--let_tautomers_change_chirality`. See `README.md` for
  important information about how Gypsum-DL treats tautomers.
-   Added Durrant-lab filters. In looking over many Gypsum-DL-generated
  variants, we have identified several substructures that, though technically
  possible, strike us as improbable. See `README.md` for examples. To discard
  molecular variants with these substructures, use the
  `--use_durrant_lab_filters` flag.
-   Rather than RDKit's PDB flavor=4, now using flavor=32.
-   PDB files now contain 2 REMARK lines describing the input SMILES string and
  the final SMILES of the ligand.
-   Added comment to `README.md` re. the need to first use drug-like filters to
  remove large molecules before Gypsum-DL processing.
-   Added comment to `README.md` re. advanced approaches for eliminating
  problematic compounds.

## [1.1.0]

-   Updated Dimorphite-DL dependency from version 1.0.0 to version 1.2.0. See
  `$PATH/gypsum_dl/gypsum_dl/Steps/SMILES/dimorphite_dl/CHANGES.md` for more
  information.
-   Updated MolVS dependency from version v0.1.0 to v0.1.1 2019 release. See
  `$PATH/gypsum_dl/gypsum_dl/molvs/CHANGELOG.md` for more information.
-   Gypsum-DL now requires mpi4py version 2.1.0 or higher. Older mpi4py versions
  can [result in deadlock if a `raise Exception` is triggered while
  multiprocessing](https://mpi4py.readthedocs.io/en/stable/mpi4py.run.html).
  Newer mpi4py versions (2.1.0 and higher) provide an alternative command line
  execution mechanism (the `-m` flag) that implements the runpy Python module.
  Gypsum-DL also requires `-m mpi4py` to run in mpi mode (e.g., `mpirun -n
  $NTASKS python -m mpi4py run_gypsum_dl.py ...-settings...`). If you experience
  deadlock, [please contact](mailto:durrantj@pitt.edu) us immediately.

  To test your version of mpi4py, open a python window and run the following
  commands:

    ```python
    >>> import mpi4py
    >>> print(mpi4py.__version__)
    3.0.1
    >>>
    ```

-   Updated the examples and documentation (`-h`) to reflect the above changes.
-   Added a Gypsum-DL citation to the print statement.

## [1.0.0]

The original version described in:

Ropp PJ, Spiegel JO, Walker JL, Green H, Morales GA, Milliken KA, Ringe JJ,
Durrant JD. Gypsum-DL: An Open-Source Program for Preparing Small-Molecule
Libraries for Structure-Based Virtual Screening. J Cheminform. 11(1):34, 2019.
[PMID: 31127411] [doi: 10.1186/s13321-019-0358-3]