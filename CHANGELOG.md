Changes
=======

1.1.5
-----

* Updated Dimorphite-DL to 1.2.4. Now better handles compounds with
  polyphosphate chains (e.g., ATP).
* Minor updates to the Durrant-lab filters:
  * When running Gypsum-DL without the `--use_durrant_lab_filters` parameter,
    Gypsum-DL now displays a warning. We strongly recommend using these
    filters, but we choose not to turn them on by default in order to maintain
    backwards compatibility.
  * Added filter to compensate for a phosphate-related bug in MolVS, one of
    Gypsum-DL's dependencies. MolVS sometimes tautomerizes `[O]P(O)([O])=O` to
    `[O][PH](=O)([O])=O`, so the Durrant-lab filters now remove any tautomers
    with substructures that match the SMARTS string `O=[PH](=O)([#8])([#8])`.
  * Added filters to compensate for frequently seen, unusual MolVS
    tautomerization of adenine. The Durrant-lab filters now remove tautomers
    with substructures that match `[#7]=C1[#7]=C[#7]C=C1` and
    `N=c1cc[#7]c[#7]1`.
  * Added filter to remove terminal iminols. While amide-iminol
    tautomerization is valid, amides are far more common, and accounting for
    this tautomerization produces many improbable iminol componds. The
    Durrant-lab filters now remove compounds with substructures that match
    `[$([NX2H1]),$([NX3H2])]=C[$([OH]),$([O-])]`.
  * Added filter to remove molecules containing `[Bi]`.
* Gypsum-DL now outputs molecules with total charges between -4e and +4e.
  Before, the cutoff was -2e to 2e. We expanded the range to permit ATP and
  other similar molecules.

1.1.4
-----

* Updated Dimorphite-DL to 1.2.3.
* Added `sys.stdout.flush()` commands to ParallelMPI.run (see
  `gypsum_dl/gypsum_dl/Parallelizer.py`) to ensure that print statements
  properly output in large MPI runs.

1.1.3
-----

* Gypsum-DL used to crash when provided with certain mal-formed SMILES
  strings. It now just skips those SMILES and warns the user that they are
  poorly formed. See Start.py:303 and MyMol.py:747.
* Durrant-lab filters now remove molecules containing metal and boron atoms.
* Some Durrant-lab filters are now applied immediately after desalting. We
  discovered that certain substructures cause Gypsum-DL to delay during the
  add-hydrogens step, specifically when Gypsum-DL generates the 3D structures
  required to rank conformers. Removing these compounds before adding
  hydrogens avoids the problem.
* Improved code formatting.
* Made minor spelling corrections to the output.

1.1.2
-----

* Bug fix: thoroughness parameter is now correctly recognized as an int when
  specified from the command line.

1.1.1
-----

* Updated Dimorphite-DL to version 1.2.2.
* Corrected spelling in user-parameter names. Parameters that previously used
  "ennumerate" now use "enumerate".
* Updated MolVS-generated tautomer filters. Previous versions of Gypsum-DL
  rejected tautomers that changed the number of _specified_ chiral centers. By
  default, Gypsum now rejects tautomers that change the total number of chiral
  centers, _both specified and unspecified_. To override the new default
  behavior (i.e., to allow tautomers that change the total number of chiral
  centers), use `--let_tautomers_change_chirality`. See `README.md` for
  important information about how Gypsum-DL treats tautomers.
* Added Durrant-lab filters. In looking over many Gypsum-DL-generated
  variants, we have identified several substructures that, though technically
  possible, strike us as improbable. See `README.md` for examples. To discard
  molecular variants with these substructures, use the
  `--use_durrant_lab_filters` flag.
* Rather than RDKit's PDB flavor=4, now using flavor=32.
* PDB files now contain 2 REMARK lines describing the input SMILES string and
  the final SMILES of the ligand.
* Added comment to `README.md` re. the need to first use drug-like filters to
  remove large molecules before Gypsum-DL processing.
* Added comment to `README.md` re. advanced approaches for eliminating
  problematic compounds.

1.1.0
-----

* Updated Dimorphite-DL dependency from version 1.0.0 to version 1.2.0. See
  `$PATH/gypsum_dl/gypsum_dl/Steps/SMILES/dimorphite_dl/CHANGES.md` for more
  information.
* Updated MolVS dependency from version v0.1.0 to v0.1.1 2019 release. See
  `$PATH/gypsum_dl/gypsum_dl/molvs/CHANGELOG.md` for more information.
* Gypsum-DL now requires mpi4py version 2.1.0 or higher. Older mpi4py versions
  can [result in deadlock if a `raise Exception` is triggered while
  multiprocessing](https://mpi4py.readthedocs.io/en/stable/mpi4py.run.html).
  Newer mpi4py versions (2.1.0 and higher) provide an alternative command line
  execution mechanism (the `-m` flag) that implements the runpy Python module.
  Gypsum-DL also requires `-m mpi4py` to run in mpi mode (e.g., `mpirun -n
  $NTASKS python -m mpi4py run_gypsum_dl.py ...-settings...`). If you
  experience deadlock, [please contact](mailto:durrantj@pitt.edu) us
  immediately.

  To test your version of mpi4py, open a python window and run the following
  commands:

    ```python
    >>> import mpi4py
    >>> print(mpi4py.__version__)
    3.0.1
    >>>
    ```

* Updated the examples and documentation (`-h`) to reflect the above changes.
* Added a Gypsum-DL citation to the print statement.

1.0.0
-----

The original version described in:

Ropp PJ, Spiegel JO, Walker JL, Green H, Morales GA, Milliken KA, Ringe JJ,
Durrant JD. Gypsum-DL: An Open-Source Program for Preparing Small-Molecule
Libraries for Structure-Based Virtual Screening. J Cheminform. 11(1):34, 2019.
[PMID: 31127411] [doi: 10.1186/s13321-019-0358-3]
