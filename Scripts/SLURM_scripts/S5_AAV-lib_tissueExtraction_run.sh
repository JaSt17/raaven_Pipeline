#!/bin/bash
#SBATCH --job-name=S5_AAV-lib_tissueExtraction_run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=01:00:00
#SBATCH -A lu2024-2-56
#SBATCH -o S5_AAV-lib_tissueExtraction_run%j.out
#SBATCH -e S5_AAV-lib_tissueExtraction_run%j.err

# exit when any command fails
set -e

# Load the required modules
module restore R_pipeline

# activate the conda environment
conda activate R_pipeline

# Run the R script
R CMD BATCH S5_AAV-lib_tissueExtraction.R S5_R_log.txt

