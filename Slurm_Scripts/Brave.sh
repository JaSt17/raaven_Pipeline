#!/bin/bash
#SBATCH --job-name=brave_analysis
#SBATCH --nodes=3
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH -A lu2024-2-8
#SBATCH -o brave_analysis%j.out
#SBATCH -e brave_analysis%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

# change the config to the BRAVE config and rerun the BRAVE pipeline
cp configs/config_BRAVE.py Python_Scripts/config.py
./Python_Scripts/S1.py
cp configs/config_BRAVE.py Python_Scripts/config.py
./Python_Scripts/S2.py
cp configs/config_BRAVE.py Python_Scripts/config.py
./Python_Scripts/S3.py
cp configs/config_BRAVE.py Python_Scripts/config.py
./Python_Scripts/S4.py
cp configs/config_BRAVE.py Python_Scripts/config.py
./Python_Scripts/S5.py
cp configs/config_BRAVE.py Python_Scripts/config.py
./Python_Scripts/S6.py
