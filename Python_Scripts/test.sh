#!/bin/bash
#SBATCH --job-name=pipeline
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=48
#SBATCH --time=12:00:00
#SBATCH -A lu2024-2-8
#SBATCH -o pipeline%j.out
#SBATCH -e pipeline%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

./Python_Scripts/S1.py
./Python_Scripts/S2.py
./Python_Scripts/S3.py
./Python_Scripts/S4.py
./Python_Scripts/S5.py
./Python_Scripts/S6.py
