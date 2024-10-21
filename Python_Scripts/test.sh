#!/bin/bash
#SBATCH --job-name=test_S2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=04:00:00
#SBATCH -A lu2024-2-56
#SBATCH -o test_S2%j.out
#SBATCH -e test_S2%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate jupyter

Python_Scripts/S2.py
