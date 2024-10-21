#!/bin/bash
#SBATCH --job-name=fastqc_multiqc
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --time=02:00:00
#SBATCH -A lu2024-2-56
#SBATCH -o qc_job_%j.out
#SBATCH -e qc_job_%j.err

# exit when any command fails
set -e

# Load the required modules
module restore fastqc_multiqc 

# Ensure output directories exist
mkdir -p ../1_quality_control/fastqc_results/
mkdir -p ../1_quality_control/multiqc_reports/

# Run FastQC on all FASTQ files
fastqc ../0_data/fastq_files/*.fastq -o ../1_quality_control/fastqc_results/ -t 48

# Run MultiQC to aggregate the FastQC results
multiqc ../1_quality_control/fastqc_results/ -o ../1_quality_control/multiqc_reports/

