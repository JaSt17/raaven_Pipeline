# Data Processing Pipeline for AAV Libraries and Tissue Samples

This Git repository contains all the code and notebooks used to process raw Illumina short-read data into a format suitable for data analysis and machine learning tasks conducted during my Master's thesis project at rAAven Therapeutics.

In general, this pipeline enables the analysis of Illumina short reads from AAV libraries and tissue samples to identify recombinant AAV variants found in different cell types.

**Note:** For better organization, each directory contains its own dedicated README file describing its contents and purpose.

---

## Table of Contents

- [Data Processing Pipeline for AAV Libraries and Tissue Samples](#data-processing-pipeline-for-aav-libraries-and-tissue-samples)
  - [Table of Contents](#table-of-contents)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
    - [Setting Up the Conda Environment](#setting-up-the-conda-environment)
    - [Installing BBmap](#installing-bbmap)
  - [Running the Pipeline](#running-the-pipeline)
    - [The Data](#the-data)
      - [Required Files](#required-files)
        - [Illumina Short-Read FASTQ Files](#illumina-short-read-fastq-files)
        - [Example-Directory](#example-directory)
        - [Additional Files](#additional-files)
        - [Example-Directory](#example-directory-1)
    - [The Config File](#the-config-file)
    - [Start the Pipeline](#start-the-pipeline)
      - [Alternative](#alternative)
      - [HPC](#hpc)
  - [Workflow of the Pipeline](#workflow-of-the-pipeline)
      - [**Step 1: Fragment Generation and Translation**](#step-1-fragment-generation-and-translation)
      - [**Step 2: Barcode and Fragment Extraction**](#step-2-barcode-and-fragment-extraction)
      - [**Step 3: Barcode Reduction and Classification**](#step-3-barcode-reduction-and-classification)
      - [**Step 4: Sample Barcode Processing**](#step-4-sample-barcode-processing)
      - [**Step 5: Data Merging**](#step-5-data-merging)
      - [**Step 6: Final Dataset Processing and Normalization**](#step-6-final-dataset-processing-and-normalization)
  - [Data Visualisation](#data-visualisation)
  - [License](#license)
  - [Acknowledgments](#acknowledgments)
  - [Contact Information](#contact-information)

---

## Prerequisites

Before running the tool, ensure the following software and libraries are installed:

- **Anaconda** (version 24.1.2)
- **BBmap** (bioinformatics tool)

---

## Installation

### Setting Up the Conda Environment

To run the pipeline, you must set up a Conda environment using the provided `pipeline.yml` file in the `envs/` directory. This can be done with the following command:

```bash
conda env create -f envs/pipeline.yml
```

**Note:** Dependencies can also be installed using a virtual environment (venv) with `pip` or manually, but Conda is recommended for consistency.

### Installing BBmap

The pipeline uses the `bbduk2` tool from the BBmap suite to extract DNA sequences from reads based on a left and a right literal. Follow these steps to install and configure BBmap:

1. Clone the BBmap repository:

```bash
   git clone https://github.com/BioInfoTools/BBMap.git
```

**or**
Download the zip file from their official website: https://sourceforge.net/projects/bbmap/

2. Navigate to the BBmap directory and ensure bbduk2 is executable:

```bash
cd BBMap
chmod +x bbduk.sh
export PATH=$PATH:$(pwd)
```

3. Verify the installation:

```bash
./bbduk2.sh --version
```

The version used in our pipeline is 39.08

**Note:** Add the BBmap directory to your PATH in your shell configuration file (e.g., .bashrc or .zshrc) to avoid repeating the export command.

---

## Running the Pipeline

### The Data

Since proprietary data from rAAven Therapeutics cannot be published, a folder named `Example/` is provided with small example files for demonstration purposes.

#### Required Files

##### Illumina Short-Read FASTQ Files

- Paired-end reads of fragments and barcodes from the AAV library.
- Barcode reads from AAV capsids to identify which variants formed capsids.
- Barcode reads from tissue samples to analyze AAV variant presence.

##### Example-Directory

- **Paired-end reads**:
  - `Example/fastq_files/R1.fastq.gz`
  - `Example/fastq_files/R2.fastq.gz`
- **AAV capsid and tissue reads**:
  - `Example/sample_fastq/AAV23_R1.fastq.gz`
  - `Example/sample_fastq/Tissue1.fastq.gz`
  - `Example/sample_fastq/Tissue2.fastq.gz`

##### Additional Files

1. **Reference Sequences**:
   - If the library was based on known protein sequences, include a reference sequence file: `reference_seq.fasta`.

2. **Codon Usage File**:
   - `wSet.csv`: Contains human-optimized codon usage percentages for *Homo sapiens*.

3. **Load List**:
   - `loadlist.csv`: Maps tissue/sample FASTQ filenames to their corresponding tissue groups.

##### Example-Directory

- **Additional Files**:
  - `Example/input/load_list.csv`
  - `Example/input/reference_seq.fasta`
  - `Example/input/wSet.csv`

### The Config File

Before running the pipeline, customize the configuration file. The config file is a Python script containing a dictionary with parameters for each pipeline step.

Config files are stored in the `configs/` subdirectory.
This directory contains all configs that were used for the different libraries in our analysis.

├── configs  
│   ├── config_BRAVE.py  
│   ├── config_Example.py  
│   ├── config_p005.py  
│   ├── config_p006.py  
│   ├── config_p007.py  
│   ├── config_undetermined_p005.py  
│   ├── config_undetermined_p006.py  
│   ├── config_undetermined_p007.py  
│   ├── README.md  


 For the example data, use the provided `configs/config_Example.py` file in this directory.

**Note:** Detailed parameter explanations are available in `configs/README.md`.

### Start the Pipeline

After setting up the configuration file and verifying that all parameters are correctly defined, you can start the pipeline. The easiest way to execute it is by running the `run_pipeline.sh` script as follows:

```bash
./Bash_scripts/run_pipeline.sh configs/example_config.py
```

In case we want to run the pipeline with an NNK library instead of one that was designed with rational design (with a reference) we can also run the code like this:

```bash
./Bash_scripts/run_pipeline.sh configs/example_config.py -nnk True
```

This way Step 2 and 3 of the pipeline will be run according of a design of a NNK library.

#### Alternative

The pipeline can of cause also be run individually by simply running the individual python scripts

```bash
./Python_Scripts/S1.py
./Python_Scripts/S2.py (S2_NNK.py)
./Python_Scripts/S3.py (S3_NNK.py)
./Python_Scripts/S4.py
./Python_Scripts/S5.py
./Python_Scripts/S6.py
```

#### HPC

In case we are working on a HPC with a Slurm queueing system we can use the Slrum scripts in the `Slurm_Scripts`directory. These scripts utulize multiple cores and this way speed up the computing time of the pipeline. Of cause the settings of the Slurm scripts can be changed accoring to the avaivlabe ressources of the HPC.

For example the complete workflow of the the Example dataset can be queued like this:

```bash
sbatch Slurm_Scripts/Example.sh
```

---

## Workflow of the Pipeline

This list provides a concise explanation of each step in the pipeline workflow. Each script processes data sequentially to generate, analyze, and organize fragments and barcodes.

---

#### **Step 1: Fragment Generation and Translation**

- Reads a FASTA file (path specified in the config file).  
- Generates all possible fragments of a specified length from nucleotide sequences.  
- Translates the reference sequence into human-optimized codon usage and saves a new reference file.  
- Creates a table with detailed fragment information and saves it as a CSV file.

---

#### **Step 2: Barcode and Fragment Extraction**

- Extracts barcodes and fragments from sequencing files.  
- Matches extracted data to the reference, pairs barcodes with fragments, and saves the results to separate files.

---

#### **Step 3: Barcode Reduction and Classification**

- Reduces barcodes using the Starcode algorithm.  
- Replaces original barcodes with Starcode-reduced versions.  
- Classifies barcodes into:
  - **Definitive barcodes**  
  - **Chimeric barcodes** (based on a threshold ratio of maximal to total read count)  
  - **Single barcodes**  
- Sets appropriate modes for each barcode type:  
  - "Single" for single-read barcodes  
  - "Def" for clean multi-read barcodes  
  - "Chimeric" for chimeric barcodes  
- Saves results in three CSV files.

---

#### **Step 4: Sample Barcode Processing**

- Extracts barcodes from provided cell samples.  
- Reduces barcodes using Starcode and matches them to corresponding fragments.  
- Processes a CSV file containing **Sample** and **Group** data.  
- Saves results in a log table and stores identified fragments in a CSV file.

---

#### **Step 5: Data Merging**

- Merges fragment information from the library with the reference fragment information.  
- Outputs a combined dataset for further analysis.

---

#### **Step 6: Final Dataset Processing and Normalization**

- Processes the identified fragments and consolidates them into a single comprehensive dataset.  
- The final dataset includes detailed information for every found fragment across different samples.  
- Normalizes read counts in the dataset to adjust for the total RNA counts in each group.  
- Outputs a normalized, ready-to-analyze dataset.

---

**Note:** A more detailed explaination of the Workflow can be found in `Supplementaries/Workflow.pdf`

## Data Visualisation

All the code required for creating data visualizations and performing analysis can be found in the `Plotting_Scripts` directory. This directory contains standalone scripts for generating specific plots, as well as Jupyter notebooks that can be executed to reproduce the visualizations step-by-step. Whether you are interested in exploring individual plots or analyzing trends across datasets, these resources provide a straightforward way to recreate and customize the visualizations.

## License

This pipeline free to use, modify, and distribute.

## Acknowledgments

I sincerely thank **rAAVen Therapeutics** for their support, guidance, and the opportunity to contribute to this project. Special thanks to my supervisor at Lund University for their mentorship and encouragement throughout this journey.

## Contact Information

If there are any problems with the pipeline please file an issue on GitHub.
