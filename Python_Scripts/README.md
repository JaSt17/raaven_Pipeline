# Python Scripts Overview

This directory contains all **Python scripts** required for the pipeline. A detailed explanation of the workflow can be found in the `Supplementaries` folder.

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

3. **S1.py, S2.py, S3.py, S4.py, S5.py, S6.py**  
   - **Purpose**  
     Contain the main code for each step (**S1** through **S6**) of the pipeline when analyzing a library created with **reference DNA sequences**.  
   - **Usage Notes**  
     - Steps must typically be run in numerical order unless specified otherwise.  
     - Each script uses the settings from `config.py` (or the specified config file in `configs/`).

4. **S2_NNK.py, S3_NNK.py, S6_NNK.py**
   - **Purpose**  
     Provide alternative versions of steps **S2** and **S3** for random **NNK libraries**, where no reference sequences exist.  

---
