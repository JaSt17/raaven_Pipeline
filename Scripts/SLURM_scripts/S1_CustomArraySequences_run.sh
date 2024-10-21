#!/bin/bash
#SBATCH --job-name=S1_CustomArraySequences_run
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=00:15:00
#SBATCH -A lu2024-2-56
#SBATCH -o S1_CustomArraySequences_run%j.out
#SBATCH -e S1_CustomArraySequences_run%j.err

# exit when any command fails
set -e

# Load the required modules
module restore R_pipeline

# Run the R script
R CMD BATCH R_Scripts/S1_CustomArraySequences.R S1_R_log.txt

