#!/bin/bash
#SBATCH --job-name=test_S5
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=48
#SBATCH --time=1:00:00
#SBATCH -A lu2024-2-8
#SBATCH -o test_S5%j.out
#SBATCH -e test_S5%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

Python_Scripts/S5.py
