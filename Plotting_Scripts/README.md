# Plotting Scripts Overview

This directory contains scripts and notebooks used for data visualization and plotting analyses.

---

## Files and Directories Description

1. **barcode_fragment_logo.ipynb**  
    - **Purpose**  
    Generates logo plots representing the barcode and fragment sequences from different libraries.  
    - **Usage Notes**  
    - Requires input sequence data.
    - Plots are saved in the `plots/logos` directory.

2. **plotting_functions.py**  
    - **Purpose**  
    Contains all custom plotting functions utilized in the `visualisation.ipynb` notebook.  

3. **visualisation.ipynb**  
   - **Purpose**  
     Main notebook for generating all visualizations (plots and tables) for data analysis.  
   - **Usage Notes**  
     - Depends on `plotting_functions.py`.
     - Saves all outputs to the `plots/` directory.

4. **plots/**  
   - **Purpose**  
     Directory for storing all generated visualization plots and output images.  
   - **Usage Notes**  
     - Automatically populated when running `barcode_fragment_logo.ipynb` or `visualisation.ipynb`.

---
