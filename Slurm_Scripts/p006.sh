#!/bin/bash
#SBATCH --job-name=p006_analysis
#SBATCH --nodes=4
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=48
#SBATCH --time=12:00:00
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
cp configs/config_p006.py Python_Scripts/config.py
# A NNK library starts at step 2 since there is no reference file
./Python_Scripts/S2_NNK.py
./Python_Scripts/S3_NNK.py
./Python_Scripts/S4.py
./Python_Scripts/S5.py
./Python_Scripts/S6_NNK.py