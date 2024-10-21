#!/bin/bash

# SRR ID passed as the first argument
srr="$1"

# Download the SRA file using prefetch
prefetch $srr

# Convert the SRA file to FASTQ format using fasterq-dump
# The --outdir option specifies the directory for FASTQ files
fasterq-dump $srr --outdir ./fastq_files --split-files

