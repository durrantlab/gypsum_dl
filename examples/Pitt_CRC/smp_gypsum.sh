#!/usr/bin/env bash
#SBATCH --job-name=gypsum_smp
#SBATCH --output=gypsum_smp_out.txt
#SBATCH --time=0:03:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cluster=smp
#SBATCH --partition=high-mem

## Define the environment
module purge
module load intel/2017.1.132 intel-mpi/2017.1.132
module load python/anaconda2.7-4.2.0

## Run the process
python run_gypsum.py -j smp_sample_molecules.json > test_smp_gypsum_output.txt

