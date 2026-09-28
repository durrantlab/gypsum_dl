#!/usr/bin/env bash
#SBATCH --job-name=gypsum_dl_mpi
#SBATCH --output=crc_mpi_gypsum_dl_out.txt
#SBATCH --time=0:07:00
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=12
#SBATCH --cluster=mpi
#SBATCH --partition=opa

## Define the environment
module purge
module load gcc/8.2.0
module load python/anaconda3.7-2018.12_westpa

## Warm the bytecode cache on a single rank before the ranks fan out. Matters
## only when running from a source checkout (an editable or path install), where
## the modules are not byte-compiled at install time and every rank would
## otherwise compile them at once on shared storage.
srun -n 1 python -c "import gypsum_dl.run"

## Run the process
mpirun -n $SLURM_NTASKS python -m mpi4py run_gypsum_dl.py -j mpi_sample_molecules.json > test_mpi_gypsum_dl_output.txt
