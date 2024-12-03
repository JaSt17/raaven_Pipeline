#!/bin/bash
#SBATCH --job-name=pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=50000
#SBATCH --time=48:00:00
#SBATCH -A lu2024-2-8
#SBATCH -o pipeline%j.out
#SBATCH -e pipeline%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

# change the config to the BRAVE config and rerun the BRAVE pipeline
cp Python_Scripts/configs/config_BRAVE.py Python_Scripts/config.py
# ./Python_Scripts/S1.py
# ./Python_Scripts/S2.py
./Python_Scripts/S3_chunked.py
./Python_Scripts/S4.py
./Python_Scripts/S5.py
./Python_Scripts/S6.py

# change the config file to p005
cp Python_Scripts/configs/config_p005.py Python_Scripts/config.py
./Python_Scripts/S1.py
./Python_Scripts/S2.py
./Python_Scripts/S3_chunked.py
./Python_Scripts/S5.py

# change the config file to p006
cp Python_Scripts/configs/config_p006.py Python_Scripts/config.py
#running only S2 and S§_NNK since p006 does not have reference data
#./Python_Scripts/S2.py
#./Python_Scripts/S3_NNK_chunked.py

# change the config file to p007
cp Python_Scripts/configs/config_p007.py Python_Scripts/config.py
./Python_Scripts/S1.py
./Python_Scripts/S2.py
./Python_Scripts/S3_chunked.py
./Python_Scripts/S5.py
