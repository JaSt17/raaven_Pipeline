#!/bin/bash
#SBATCH --job-name=S4_Fragment_run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=12:00:00
#SBATCH -A lu2024-2-56
#SBATCH -o S4_Fragment_run%j.out
#SBATCH -e S4_Fragment_run%j.err

# exit when any command fails
set -e

# Load the required modules
module restore R_pipeline

# activate the conda environment
conda activate R_pipeline

# Run the R script
R CMD BATCH S4_Fragment_translation.R S4_R_log.txt

