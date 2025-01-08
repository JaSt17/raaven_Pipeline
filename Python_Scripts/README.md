# Python Scripts Overview

This directory contains all **Python scripts** required for the pipeline. A detailed explanation of the workflow can be found in both PDF and TXT formats in the `Supplementaries` folder.

---

## Scripts Description

1. **config.py**  
   - **Purpose**  
     Holds all configuration parameters needed to run the pipeline.  
   - **Usage Notes**  
     - The actual configuration files are stored in the `configs/` directory.  
     - **Important**: The config must be updated before running any individual step of the pipeline.  
     - If you use the Slurm scripts or `Bash_Scripts/run_pipeline.sh`, the configuration is automatically adjusted.

2. **costum_functions.py**  
   - **Purpose**  
     A collection of helper functions that are needed in multiple steps of the pipeline.  
   - **Why It’s Helpful**  
     - Prevents code duplication by allowing you to import common functions rather than rewriting them.  

3. **S1.py, S2.py, S3.py, S4.py, S5.py, S6.py**  
   - **Purpose**  
     Contain the main code for each step (**S1** through **S6**) of the pipeline when analyzing a library created with **reference DNA sequences**.  
   - **Usage Notes**  
     - Steps must typically be run in numerical order unless specified otherwise.  
     - Each script uses the settings from `config.py` (or the specified config file in `configs/`).

4. **S2_NNK.py, S3_NNK.py**  
   - **Purpose**  
     Provide alternative versions of steps **S2** and **S3** for random **NNK libraries**, where no reference sequences exist.  
   - **Why It’s Helpful**  
     - Offers flexibility for **randomized** library analysis.  
     - Streamlines workflows when you don’t have a predefined set of reference sequences.

---

> **Note**: For more in-depth explanations, please see the **PDF** and **TXT** documentation in the `Supplementaries` folder.
