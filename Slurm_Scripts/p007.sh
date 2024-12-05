#!/bin/bash
#SBATCH --job-name=p007_analysis
#SBATCH --nodes=5
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=48
#SBATCH --time=48:00:00
#SBATCH -A lu2024-2-8
#SBATCH -o p007_analysis%j.out
#SBATCH -e p007_analysis%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

# change the config file to p007
cp Python_Scripts/configs/config_p007.py Python_Scripts/config.py
./Python_Scripts/S1.py
./Python_Scripts/S2.py
./Python_Scripts/S3.py
./Python_Scripts/S5.py