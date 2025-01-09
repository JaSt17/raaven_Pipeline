#!/bin/bash
#SBATCH --job-name=example_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=1:00:00
#SBATCH -A lu2024-2-8
#SBATCH -o example_analysis%j.out
#SBATCH -e example_analysis%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

# change the config to the Example config and rerun the pipeline
cp configs/config_Example.py Python_Scripts/config.py
./Python_Scripts/S1.py
./Python_Scripts/S2.py
./Python_Scripts/S3.py
./Python_Scripts/S4.py
./Python_Scripts/S5.py
./Python_Scripts/S6.py
