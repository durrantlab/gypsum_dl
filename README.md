# Gypsum-DL 1.2.1

Gypsum-DL is a free, open-source program for preparing 3D small-molecule
models. Beyond simply assigning atomic coordinates, Gypsum-DL accounts for
alternate ionization, tautomeric, chiral, cis/trans isomeric, and
ring-conformational forms. It is released under the Apache License, Version
2.0 (see `LICENSE.txt`).

## Citation

If you use Gypsum-DL in your research, please cite:

Ropp, Patrick J., Jacob O. Spiegel, Jennifer L. Walker, Harrison Green,
Guillermo A. Morales, Katherine A. Milliken, John J. Ringe, and Jacob D.
Durrant. (2019) "Gypsum-DL: An Open-source Program for Preparing
Small-molecule Libraries for Structure-based Virtual Screening." Journal of
Cheminformatics 11:1. doi:10.1186/s13321-019-0358-3.

Ropp PJ, Kaminsky JC, Yablonski S, Durrant JD (2019) Dimorphite-DL: An
open-source program for enumerating the ionization states of drug-like small
molecules. J Cheminform 11:14. doi:10.1186/s13321-019-0336-9.

## Getting Started

To run Gypsum-DL, acquire a copy of this repository, either by git clone or by
download. Install the required dependencies via your favorite python package
manager. We suggest using Anaconda to manage packages:

```bash
conda install -c rdkit rdkit numpy scipy mpi4py
```

If you encounter an error like `ImportError: libboost_python3.so.1.65.1:
cannot open shared object file: No such file or directory`, [a helpful
user](https://durrantlab.pitt.edu/forums/topic/importerror-libboost_python3-so-1-65-1-solved/)
suggested installing RDKit this way instead:

```bash
conda create -c conda-forge --name gypsum_dl_env rdkit numpy scipy mpi4py -y
conda activate gypsum_dl_env
```

## Command-Line Parameters

Gypsum-DL accepts the following command-line parameters:

```text
  -h, --help            show this help message and exit
  --json param.json, -j param.json
                        Name of a json file containing all parameters.
                        Overrides all other arguments specified at the
                        commandline.
  --source input.smi, -s input.smi
                        Name of the source file (e.g., input.smi). Note:
                        support for SMI (SMILES) files is better than support
                        for SDF files, though Gypsum-DL can handle both.
  --output_folder OUTPUT_FOLDER, -o OUTPUT_FOLDER
                        The path to an existing folder where the Gypsum-DL
                        output file(s) will be saved.
  --job_manager {mpi,multiprocessing,serial}
                        Determine what style of multiprocessing to use: mpi,
                        multiprocessing, or serial. Serial will override the
                        num_processors flag, forcing it to be one. MPI mode
                        requires mpi4py 2.1.0 or higher and should be executed
                        as: mpirun -n $NTASKS python -m mpi4py
                        run_gypsum_dl.py ...-settings...
  --num_processors N, -p N
                        Number of processors to use for parallel calculations.
  --max_variants_per_compound V, -m V
                        The maximum number of variants to create per input
                        molecule.
  --thoroughness THOROUGHNESS, -t THOROUGHNESS
                        How widely to search for low-energy conformers. Larger
                        values increase run times but can produce better
                        results.
  --separate_output_files
                        Indicates that the outputs should be split between
                        files. If true, each output .sdf file will correspond
                        to a single input file, but different 3D conformers
                        will still be stored in the same file.
  --add_pdb_output      Indicates that the outputs should also be written in
                        the .pdb format. Creates one PDB file for each
                        molecular variant.
  --add_html_output     Indicates that the outputs should also be written in
                        the .html format, for debugging. Attempts to open a
                        browser for viewing.
  --min_ph MIN          Minimum pH to consider.
  --max_ph MAX          Maximum pH to consider.
  --pka_precision D     Size of pH substructure ranges. See Dimorphite-DL
                        publication for details.
  --skip_optimize_geometry
                        Skips the optimization step.
  --skip_alternate_ring_conformations
                        Skips the non-aromatic ring-conformation generation
                        step.
  --skip_adding_hydrogen
                        Skips the ionization step.
  --skip_making_tautomers
                        Skips tautomer-generation step.
  --skip_enumerate_chiral_mol
                        Skips the ennumeration of unspecified chiral centers.
  --skip_enumerate_double_bonds
                        Skips the ennumeration of double bonds.
  --let_tautomers_change_chirality
                        Allow tautomers that change the total number of chiral
                        centers (see README.md for further explanation).
  --use_durrant_lab_filters
                        Use substructure filters to remove molecular variants
                        that, though technically possible, were judged
                        improbable by members of the Durrant lab. See
                        README.md for more details.
  --2d_output_only      Skips the generate-3D-models step.
  --cache_prerun, -c    Run this before running Gypsum-DL in mpi mode.
  --test                Tests Gypsum-DL to check for programming bugs.
```

## Examples of Use

Prepare a virtual library and save all 3D models to a single SDF file in the
present directory:

```bash
python run_gypsum_dl.py --source ./examples/sample_molecules.smi
```

Instead save all 3D models to a different, existing folder:

```bash
python run_gypsum_dl.py --source ./examples/sample_molecules.smi \
   --output_folder /my/folder/
```

Additionally save the models associated with each input molecule to separate
files:

```bash
python run_gypsum_dl.py --source ./examples/sample_molecules.smi \
    --output_folder /my/folder/ --separate_output_files
```

In addition to saving a 3D SDF file, also save 3D PDB files and an HTML file
with 2D structures (for debugging).

```bash
python run_gypsum_dl.py --source ./examples/sample_molecules.smi \
    --output_folder /my/folder/ --add_pdb_output --add_html_output
```

Save at most two variants per input molecule:

```bash
python run_gypsum_dl.py --source ./examples/sample_molecules.smi \
    --output_folder /my/folder/ --max_variants_per_compound 2
```

Control how Gypsum-DL ionizes the input molecules:

```bash
python run_gypsum_dl.py --source ./examples/sample_molecules.smi \
    --output_folder /my/folder/ --min_ph 12 --max_ph 14 --pka_precision 1
```

Run Gypsum-DL in serial mode (using only one processor):

```bash
python run_gypsum_dl.py --source ./examples/sample_molecules.smi \
    --job_manager serial
```

Run Gypsum-DL in multiprocessing mode, using 4 processors:

```bash
python run_gypsum_dl.py --source ./examples/sample_molecules.smi \
    --job_manager multiprocessing --num_processors 4
```

Run Gypsum-DL in mpi mode using all available processors:

```bash
mpirun -n $NTASKS python -m mpi4py  run_gypsum_dl.py --source ./examples/sample_molecules.smi \
    --job_manager mpi --num_processors -1
```

Gypsum-DL can also take parameters from a JSON file:

```bash
python run_gypsum_dl.py --json myparams.json
```

Where `myparams.json` might look like:

```json
{
    "source": "./examples/sample_molecules.smi",
    "separate_output_files": true,
    "job_manager": "multiprocessing",
    "output_folder": "/my/folder/",
    "add_pdb_output": true,
    "add_html_output": true,
    "num_processors": -1
}
```

## Important Caveats

### Large Molecules

Gypsum-DL is designed to process drug-like molecules. Generating 3D structures
for larger molecules takes a very long time. For example, in our tests it
takes Gypsum-DL a very long time to process this molecule:
`CCCC[C@@H](C(N[C@H]1CC(NCCCC[C@H](NC([C@@H](NC([C@@H](NC([C@@H](NC([C@@H]2CCCN2C1=O)=O)Cc3ccccc3)=O)CCCNC(N)=N)=O)Cc(c[nH]4)c5c4cccc5)=O)C(N6CCC[C@H]6C(N[C@@H](C(C)C)C(N)=O)=O)=O)=O)=O)NC(C)=O`

You may wish to run your compounds through a drug-like filter before
processing them with Gypsum-DL.

### Tautomers

Gypsum-DL uses MolVS to generate tautomers. While MolVS is effective, we have
noticed that it sometimes generates inappropriate tautomers that change the
total number of chiral centers, e.g. `O=C(c1ccc(CN)cc1)N` to
`N=Cc1ccc(C(O)N)cc1`. But some legitimate tautomers also change the number of
chiral centers, e.g., `C[C@@H](C(C)=O)F` to `C/C(F)=C(C)\O`.

To compensate for this MolVS bug, by default Gypsum-DL rejects all tautomers
that change the total number of chiral centers. Use the
`--let_tautomers_change_chirality` flag if you would like to retain these
tautomers instead. As always, be sure to examine the structures that Gypsum-DL
outputs to ensure they are chemically feasible.

### Durrant-Lab Filters

In looking over many Gypsum-DL-generated variants, we have identified a number
of substructures that, though technically possible, strike us as improbable or
otherwise poorly suited for virtual screening. Here are some examples:

* `C=[N-]`
* `[N-]C=[N+]`
* `[nH+]c[n-]`
* `[#7+]~[#7+]`
* `[#7-]~[#7-]`
* `[!#7]~[#7+]~[#7-]~[!#7]`
* `[#5]` (boron)
* `O=[PH](=O)([#8])([#8])`
* `N=c1cc[#7]c[#7]1`
* `[$([NX2H1]),$([NX3H2])]=C[$([OH]),$([O-])]`
* Metals

If you'd like to discard molecular variants with substructures such as these,
use the `--use_durrant_lab_filters` flag.

### Highly Constrained Ring Systems

Some users have reported that Gypsum-DL fails to produce 3D models when
processing molecules with highly constrained ring systems, such as amantadine
compounds (e.g., this molecule from ChemBridge:
`CC1=CC=CN2N=CC(C(=O)NC34CC5CC(C3)CC(C5)(C4)N3C=NC=N3)=C12`). Increasing the
`thoroughness` parameter may help in these cases.

### Memory Considerations

When processing large libraries, Gypsum-DL requires substantial memory. Some
users have reported that the program suddenly stops in these situations. To
correct the problem, either increase the available memory, or divide your
library into several smaller files and processes them sequentially.

### Advanced Methods for Eliminating Problematic Compounds

Gypsum-DL aims to enumerate many possible variant forms, including forms that
are not necessarily probable. Beyond applying Durrant-Lab filters, several
methods allow users to exclude other potentially problematic forms:

1. Identify the steps Gypsum-DL takes to generate a given problematic form
   (see the "Genealogy" field of every output SDF file). Then use parameters
   such as `--skip_optimize_geometry`, `--skip_alternate_ring_conformations`,
   `--skip_adding_hydrogen`, `--skip_making_tautomers`,
   `--skip_enumerate_chiral_mol`, or `--skip_enumerate_double_bonds` to skip
   the problem-causing step. This fix is easy, but it may unexpectedly impact
   unrelated compounds.
2. Consider adjusting the `--min_ph`, `--max_ph`, or `--pka_precision`
   parameters if Gypsum-DL is producing compounds with undesired protonation
   states. Alternatively, you can delete specific protonation rules by
   modifying the
   `gypsum_dl/Steps/SMILES/dimorphite_dl/site_substructures.smarts` file.
3. Add to the Durrant-Lab filters if there is a specific substructure you
   would like to avoid (e.g., imidic acid due to amide/imidic-acid
   tautomerization). Simplify modify the
   `gypsum_dl/Steps/SMILES/DurrantLabFilter.py` file.
