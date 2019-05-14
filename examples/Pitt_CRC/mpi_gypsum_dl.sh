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
module load intel/2017.1.132 intel-mpi/2017.1.132
module load python/anaconda2.7-4.2.0

## Run the process
mpirun -n 1 python -m mpi4py run_gypsum_dl.py -c
mpirun -n $SLURM_NTASKS python -m mpi4py run_gypsum_dl.py -j mpi_sample_molecules.json > test_mpi_gypsum_dl_output.txt
