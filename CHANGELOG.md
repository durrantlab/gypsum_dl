Changes
=======

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
