#!/bin/bash
#SBATCH --job-name=brave_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
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
# we are using the NNK S2 script for BRAVE since we are looking for multiple different fragment types (14aa, 22aa, etc)
./Python_Scripts/S2_NNK.py
./Python_Scripts/S3.py
./Python_Scripts/S4.py
./Python_Scripts/S5.py
./Python_Scripts/S6.py
