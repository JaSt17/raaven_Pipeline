#!/bin/bash
#SBATCH --job-name=S2_AAV-lib_run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:15:00
#SBATCH -A lu2024-2-56
#SBATCH -o S2_AAV-lib_run%j.out
#SBATCH -e S2_AAV-lib_run%j.err

# exit when any command fails
set -e

# Load the required modules
module restore R_pipeline

# Run the R script
R CMD BATCH S2_AAV-lib_libraryExtraction.R S2_R_log.txt

