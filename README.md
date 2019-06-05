# Gypsum-DL 1.1.0

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

## Getting Started

To run Gypsum-DL, acquire a copy of this repository, either by git clone or by
download. Install the required dependencies via your favorite python package
manager. We suggest using Anaconda to manage packages:

```bash
conda install -c rdkit rdkit numpy scipy mpi4py
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
                        Name of the source file (e.g., input.smi).
  --output_folder OUTPUT_FOLDER, -o OUTPUT_FOLDER
                        The path to an existing folder where the Gypsum-DL
                        output file(s) will be saved.
  --job_manager {mpi,multiprocessing,serial}
                        Determine what style of multiprocessing to use: mpi,
                        multiprocessing, or serial. Serial will override the
                        num_processors flag, forcing it to be one. MPI mode
                        requires mpi4py 2.1.0 or higher and should be executed
                        as:
                            mpirun -n $NTASKS python -m mpi4py run_gypsum_dl.py \
                            ...-settings...
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
  --skip_ennumerate_chiral_mol
                        Skips the ennumeration of unspecified chiral centers.
  --skip_ennumerate_double_bonds
                        Skips the ennumeration of double bonds.
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
