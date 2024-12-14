#!/bin/bash
#SBATCH --job-name=all_analysis
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=96:00:00
#SBATCH -A lu2024-2-8
#SBATCH -o all_analysis%j.out
#SBATCH -e all_analysis%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

# change the config to the BRAVE config and rerun the BRAVE pipeline
cp Python_Scripts/configs/config_BRAVE.py Python_Scripts/config.py
./Python_Scripts/S1.py
./Python_Scripts/S2.py
./Python_Scripts/S3.py
./Python_Scripts/S4.py
./Python_Scripts/S5.py
./Python_Scripts/S6.py

# change the config file to p005
cp Python_Scripts/configs/config_p005.py Python_Scripts/config.py
#./Python_Scripts/S1.py
#./Python_Scripts/S2.py
./Python_Scripts/S3.py
./Python_Scripts/S4.py
./Python_Scripts/S5.py

# change the config file to p006
cp Python_Scripts/configs/config_p006.py Python_Scripts/config.py
#running only S2 and S3_NNK since p006 does not have reference data
#./Python_Scripts/S2.py
./Python_Scripts/S3_NNK.py
./Python_Scripts/S4.py

# change the config file to p007
cp Python_Scripts/configs/config_p007.py Python_Scripts/config.py
#./Python_Scripts/S1.py
#./Python_Scripts/S2.py
./Python_Scripts/S3.py
./Python_Scripts/S4.py
./Python_Scripts/S5.py