#!/usr/bin/env bash
#SBATCH --job-name=gypsum_mpi
#SBATCH --output=crc_mpi_gypsum_out.txt
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
srun -n $SLURM_NTASKS --mpi=pmi2 python run_gypsum.py -j mpi_sample_molecules.json > test_mpi_gypsum_output.txt

