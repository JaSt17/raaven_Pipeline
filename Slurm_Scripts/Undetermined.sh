#!/bin/bash
#SBATCH --job-name=p006_analysis
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

