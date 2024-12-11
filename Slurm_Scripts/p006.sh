#!/bin/bash
#SBATCH --job-name=p006_analysis
#SBATCH --nodes=12
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=48
#SBATCH --time=8:00:00
#SBATCH -A lu2024-2-8
#SBATCH -o p006_analysis%j.out
#SBATCH -e p006_analysis%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

# change the config file to p006
cp Python_Scripts/configs/config_p006.py Python_Scripts/config.py
#running only S2 and S3_NNK since p006 does not have reference data
#./Python_Scripts/S2.py
./Python_Scripts/S3_NNK.py