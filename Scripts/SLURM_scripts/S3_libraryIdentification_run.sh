#!/bin/bash
#SBATCH --job-name=S3_libraryIdentification_run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:15:00
#SBATCH -A lu2024-2-56
#SBATCH -o S3_libraryIdentification_run%j.out
#SBATCH -e S3_libraryIdentification_run%j.err

# exit when any command fails
set -e

# Load the required modules
module restore R_pipeline

# activate the conda environment
conda activate R_pipeline

# Run the R script
R CMD BATCH S3_libraryIdentification.R S3_R_log.txt

