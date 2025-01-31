#!/bin/bash
#SBATCH --job-name=undetermined_analysis
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=48
#SBATCH --time=24:00:00
#SBATCH -A lu2024-2-8
#SBATCH -o undetermined_analysis%j.out
#SBATCH -e undetermined_analysis%j.err

# exit when any command fails
set -e

module load Anaconda3/2024.02-1

# activate the conda environment
source config_conda.sh
conda activate python_pipeline

# for p005
# change the config filte to undetermined_p005.py
cp configs/config_undetermined_p005.py Python_Scripts/config.py
# run the second step of the pipeline
./Python_Scripts/S2.py
# rename the S2.log file to S2_undetermined_p005.log
mv raav-60/undetermined/logs/S2.log raav-60/undetermined/logs/S2_undetermined_p005.log

# for p006
# change the config filte to undetermined_p006.py
cp configs/config_undetermined_p006.py Python_Scripts/config.py
# run the second step of the pipeline
./Python_Scripts/S2_NNK.py
# rename the S2.log file to S2_undetermined_p006.log
mv raav-60/undetermined/logs/S2.log raav-60/undetermined/logs/S2_undetermined_p006.log

# for p007
# change the config filte to undetermined_p007.py
cp configs/config_undetermined_p007.py Python_Scripts/config.py
# run the second step of the pipeline
./Python_Scripts/S2.py
# rename the S2.log file to S2_undetermined_p007.log
mv raav-60/undetermined/logs/S2.log raav-60/undetermined/logs/S2_undetermined_p007.log

