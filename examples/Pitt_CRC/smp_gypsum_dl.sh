#!/usr/bin/env bash
#SBATCH --job-name=gypsum_dl_smp
#SBATCH --output=gypsum_dl_smp_out.txt
#SBATCH --time=0:03:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cluster=smp
#SBATCH --partition=high-mem

## Define the environment
module purge
module load gcc/8.2.0
module load python/anaconda3.7-2018.12_westpa

## Run the process
python run_gypsum_dl.py -j smp_sample_molecules.json > test_smp_gypsum_dl_output.txt
