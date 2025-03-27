# Plotting Scripts Overview

This directory contains scripts and Jupyter notebooks used for generating data visualizations and plots related to sequence analysis and motif exploration across different tissue types.

---

## Files and Directories Description

### 1. **motifs/**
- **Purpose**  
  Stores the output motif files generated using the [MEME Suite](https://meme-suite.org/meme/), a widely-used suite of tools for discovering statistically significant sequence motifs in sets of biological sequences.

- **Details**  
  We used MEME as an exploratory tool to identify recurring patterns or motifs in amino acid sequences that are associated with specific tissue-targeting fragments. By analyzing libraries grouped by tissue type, MEME helped us uncover conserved regions that might be functionally relevant or involved in tissue specificity.

- **Access Note**  
  Due to confidentiality reasons, the contents of the `motifs/` directory are not publicly available.
---

### 2. **barcode_fragment_logo.ipynb**  
- **Purpose**  
  Generates sequence logo plots for barcode and fragment sequences from different libraries.

- **Usage Notes**  
  - Requires input sequence data.
  - Plots are saved in the `plots/logos` directory.

---

### 3. **plotting_functions.py**  
- **Purpose**  
  Contains custom plotting functions used across notebooks, especially within `visualisation.ipynb`.

---

### 4. **visualisation.ipynb**  
- **Purpose**  
  Main notebook for generating visual summaries (plots and tables) for data analysis.

- **Usage Notes**  
  - Depends on `plotting_functions.py`.
  - Outputs are saved to the `plots/` directory.

---

### 5. **plots/**  
- **Purpose**  
  Directory for storing all generated plots and visual outputs.

- **Usage Notes**  
  - Automatically populated when running `barcode_fragment_logo.ipynb` or `visualisation.ipynb`.

---
