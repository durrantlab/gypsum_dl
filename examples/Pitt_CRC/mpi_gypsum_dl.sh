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
srun -n 1 --mpi=pmi2 python run_gypsum_dl.py -c
srun -n $SLURM_NTASKS --mpi=pmi2 python run_gypsum_dl.py -j mpi_sample_molecules.json > test_mpi_gypsum_dl_output.txt