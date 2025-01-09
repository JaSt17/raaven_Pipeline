# Slurm Scripts Overview

This directory contains **Slurm scripts** that can be used to run the **pipeline** for different libraries. Each script holds the workflow for a specific library (or set of libraries). Below is a brief explanation of the structure and purpose of these scripts.

---

## Slurm Header

The Slurm header sets the resources that will be used on the HPC and provides basic job information:

- **--job-name**: Sets the name displayed in the Slurm queue.
- **--nodes/--ntasks/--cpus-per-task**: Define the resources (nodes, tasks, CPUs) allocated for the job.
- **--time**: Specifies the maximum wall time for the job. If the job exceeds this time, it is terminated.
- **-A**: Specifies the account to be charged for computing hours.
- **-o / -e**: Direct standard output/error (`stdout`/`stderr`) to files. `%j` is replaced by the job ID in the file names.

---

## Preparation

Each Slurm script first loads the **Anaconda3** module and activates the **conda environment** needed to run the various steps of the pipeline.

---

## Scripts Description

This directory contains six Slurm scripts, each named after the dataset it analyzes:

- **Brave**
- **Example**
- **p005**
- **p006**
- **p007**
- **Undetermined**

The first five scripts (**Brave**, **Example**, **p005**, **p006**, **p007**) can all be used to analyze their respective datasets. 

The **Undetermined** script runs only the second step of the pipeline on any undetermined reads. These are the barcodes and fragments from the **p005**, **p006**, and **p007** pipelines that could not be distinguished by the sequencer.
