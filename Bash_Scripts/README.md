# Bash Scripts Overview

This directory contains **Bash scripts** that were used for  **intermediate steps** in the data analysis of different libraries. Each script has a specific purpose to streamline the workflow or combine data from multiple sources.

---

## Script Descriptions

1. **combine_undetermined_with_libraries.sh**  
   - **Purpose**  
     Concatenates barcodes and fragments from the three libraries (**p005**, **p006**, and **p007**) with those found in the undetermined sequencing reads.  
   - **Why It’s Important**  
     - Ensures that undetermined reads are integrated into the analysis.  
     - Helps achieve a more **complete** understanding of the data by not discarding potentially valuable information.  
   - **Usage Notes**  
     - Run this script **after** you’ve obtained the barcodes/fragments for each library and the undetermined reads.  
     - Produces a combined dataset that can be used for subsequent analysis steps.

2. **run_pipeline.sh**  
   - **Purpose**  
     Executes the **entire pipeline** for a given configuration script in one go, rather than running each step manually.  
   - **Usage Notes**  
     - Best suited for **smaller datasets** or quick tests.  
     - For **large datasets**, consider using **High Performance Computing (HPC)** with **Slurm** scripts to efficiently utilize computing resources and avoid potential bottlenecks.

---

> **Note**: If you are working with large volumes of data, it is recommended to transition to an HPC environment or utilize Slurm for reliable and scalable performance.
