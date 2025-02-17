#!/bin/bash
#SBATCH --job-name=brave_22aa_analysis
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH -A lu2024-2-79
#SBATCH -o %x_%j.out  # %x = job name, %j = job ID
#SBATCH -e %x_%j.err  # %x = job name, %j = job ID

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

# change the config to the BRAVE config and rerun the BRAVE pipeline
cp configs/config_BRAVE_22aa.py Python_Scripts/config.py
#./Python_Scripts/S1.py
cp configs/config_BRAVE_22aa.py Python_Scripts/config.py
#./Python_Scripts/S2.py
cp configs/config_BRAVE_22aa.py Python_Scripts/config.py
./Python_Scripts/S3.py
cp configs/config_BRAVE_22aa.py Python_Scripts/config.py
./Python_Scripts/S4.py
cp configs/config_BRAVE_22aa.py Python_Scripts/config.py
./Python_Scripts/S5.py
cp configs/config_BRAVE_22aa.py Python_Scripts/config.py
./Python_Scripts/S6.py
