# Plotting Scripts Overview

This directory contains scripts and Jupyter notebooks used for generating data visualizations and plots related to sequence analysis and motif exploration across different tissue types.

---

## Files and Directories Description

### 1. **3D_protein_plot**
- **Purpose**  
  Contains the scripts and visualization files used to generate 3D structural illustrations of the VP1 protein and the full AAV capsid. These visualizations are intended to support the interpretation of sequence-based findings by placing them in a structural context.

- **Details**  
  This module uses **py3Dmol**, a Python wrapper for interactive 3D molecular visualization in the browser, to highlight the peptide insertion site on the AAV capsid structure. Specifically, it focuses on the VP1 subunit to illustrate the spatial location of experimentally inserted motifs.

---

### 2. **Library_Analysis**
- **Purpose**  
  Contains the scripts used to analyze the composition and overlap of peptide fragments across different libraries.

- **Details**  
  This directory includes tools for identifying identical between libraries, helping assess library diversity and cross-contamination. It also supports the comparison of multiple extraction and clustering strategies, such as using **Starcode** with or without different parameters to control for **chimeric barcodes**.

---

### 3. **barcode_fragment_logo.ipynb**  
- **Purpose**  
  Generates sequence logo plots for barcode and fragment sequences from different libraries.

- **Usage Notes**  
  - Requires input sequence data.
  - Plots are saved in the `plots/logos` directory.

---

### 4. **plotting_functions.py**  
- **Purpose**  
  Contains custom plotting functions used across notebooks, especially within `visualisation.ipynb`.

---

### 5. **visualisation.ipynb**  
- **Purpose**  
  Main notebook for generating visual summaries (plots and tables) for data analysis.

- **Usage Notes**  
  - Depends on `plotting_functions.py`.
  - Outputs are saved to the `plots/` directory.

---

### 6. **plots/**  
- **Purpose**  
  Directory for storing all generated plots and visual outputs.

- **Usage Notes**  
  - Automatically populated when running `barcode_fragment_logo.ipynb` or `visualisation.ipynb`.

---
