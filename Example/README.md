# Example Data to Explore the Workflow

Since the complete sequencing data from rAAVen Therapeutics is confidential, we have provided a **small example dataset** to explore the workflow and understand how it operates.

## Config File for the Example

- The config file used to run the pipeline on this example data is stored in:  
  `configs/config_Example`.

This README gives a brief explanation of the different files and directories in this example output.

---

## Directory & File Overview

1. **barcode_fragments/**  
   - Contains the paired fragments and barcodes found from the library sequencing data.

2. **fastq_files/**  
   - Contains the **raw sequencing data** of the library that should be analyzed.

3. **found_barcodes/**  
   - Contains a CSV file for every sample/tissue that was searched for barcodes.
   - Each file lists all barcodes with their corresponding fragments from the library and includes all necessary fragment information.

4. **input/**  
   - Holds all input data necessary to run the pipeline.
   - The required data is explained in the main project’s README file.

5. **logs/**  
   - Contains all **log files** from every step of the pipeline.
   - These logs provide detailed information about each processing step.

6. **sample_fastq/**  
   - Holds all **sample files** that should be analyzed for barcodes.

7. **final_fragments_summary.csv**  
   - Final output of the pipeline. Shows all found fragments in all tissues as well as in the library.

8. **found_barcodes_report.csv**  
   - A brief summary showing the **number of barcodes** found in each sample.
   - Also indicates how many barcodes remained after barcode reduction.

9. **library_barcodes.csv**  
   - Three files containing the **output of S3** of the pipeline:
     1. **Single barcodes** (appeared only once)  
     2. **Chimeric barcodes** (matched to multiple fragments)  
     3. **Valid barcodes** (used in subsequent pipeline steps)

10. **LUT.csv**  
    - The Lookup Table (LUT) containing information about **every single fragment** created from the reference file.

11. **pos_library_barcodes/**  
    - Provides information about each fragment created from the reference sequence that was found in the library.
    - Includes all barcodes associated with each fragment and corresponding counts.

12. **SortedFragments.txt**  
    - A text file listing **every single fragment** created in **S1** of the pipeline.

---
