# 📘 User Guide: How to Use the Analysis Pipeline

Welcome to the documentation for the analysis pipeline. This guide will walk you through how to run the pipeline, analyze new sequencing data, and manage your mRNA sample groups. Examples and common use cases are included to help you get started quickly.

## How To Sections

- [How to Run the Pipeline](#how-to-run-the-pipeline)
- [How to Run the Analysis for a New Library (e.g., New PacBio Run)](#how-to-run-the-analysis-for-a-new-library-eg-new-pacbio-run)
- [How to Add New mRNA Samples for Analysis](#how-to-add-new-mrna-samples-for-analysis)
- [How to Add a New Group to the mRNA Analysis](#how-to-add-a-new-group-to-the-mrna-analysis)
- [How to Open the UI with Streamlit](#how-to-open-the-ui-with-streamlit)
- [Troubleshooting & Tips](#-troubleshooting--tips)

---

## How to Run the Pipeline

The pipeline can be executed in multiple ways depending on your compute environment:

### Option 1: Using SLURM (on a cluster)

Each project contains a `Slurm_scripts` folder. Inside, you will find SLURM submission scripts to run the pipeline:

```bash
sbatch Slurm_scripts/library.sh
```

### Option 2: Direct Bash Execution (e.g., on the Monster machine)

You can also run the SLURM script like a regular bash script:

```bash
bash Slurm_scripts/library.sh
```

### Option 3: Run Individual Pipeline Steps

To run specific steps:

1. Copy the corresponding config file for the step you want to run into the `Python_scripts` folder.
2. Then execute the script manually. For example:

```bash
cp Projects/Bluejay/Libraries/p034/Plasmid_run_1/config.py Python_Scripts/config.py
./Python_Scripts/S2.py
```

---

## How to Run the Analysis for a New Library (e.g., New PacBio Run)

Follow these steps to process a new sequencing run:

### 1. Prepare the Raw Sequencing Data

- Navigate to the `Seq_Data/` folder of your project (e.g., `Bluejay/Seq_Data/`).
- Create a new folder for the new sequencing run, e.g.:

```bash
mkdir Bluejay/Seq_Data/PacBio_run_2
```

- Place the FASTQ file(s) inside this folder:
  - For **PacBio (single-end)**: use `library.fastq.gz` (e.g. `p034.fastq.gz`)
  - For **Illumina (paired-end)**: use `library_R1.fastq.gz` and `library_R2.fastq.gz` (e.g. `p034_R1.fastq.gz` & `p034_R2.fastq.gz`)

### 2. Link to the Library Folder

Navigate to the relevant library directory:

```bash
cd Library/p034
```

- Create a new folder here **with the exact same name** as the sequencing folder created above:

```bash
mkdir Library/p034/PacBio_run_2
```

### 3. Generate or Update the Config File

Run the following script:

```bash
cd Project/Bluejay/Libraries/
./update_configs.py
```

This should create a new config file for the run.

**Note:** If this is the first time processing this library, you'll need to **manually edit the config file** to set correct literals for fragment/barcode extraction.

### 4. Run the Pipeline

You can now:

- Copy the new config to the `Python_scripts` folder and run each step manually.
- OR add it to a SLURM script and submit it to the cluster.

---

## How to Add New mRNA Samples for Analysis

### 1. Add the Sample File

Place the new mRNA sample FASTQ file(s) in the `Seq_Data/Samples/` folder.

### 2. Update the Annotation CSV

Depending on the fluorophore:

- For GFP samples, edit `annotation_GFP.csv`
- For mCherry samples, edit `annotation_mCherry.csv`

Add a new line with:
| Filename | Group         |
|------------|------------------|
| `GFP_A01_S1_R1_combined.fastq.gz`|`7601_CD_1`      |


### 3. Add to a Group (Optional)

Edit the appropriate group annotation file:

- `annotation_Groups_GFP.csv` or `annotation_Groups_mCherry.csv`

Include the sample in an existing group or define a new one (see next section for details).

### 4. Re-run the Analysis

After updating the annotations you only have to rerun step 4 to 6 of the pipeline:

```bash
./Python_scripts/S4.py
./Python_scripts/S5.py
./Python_scripts/S6.py
```

Your new sample will be automatically included in the pipeline.

---

## How to Add a New Group to the mRNA Analysis

To group mRNA samples for comparative analysis:

### 1. Open the Group Annotation File

Navigate to:

```bash
Seq_Data/Samples/annotation_Groups_GFP.csv
```
or 
```bash
Seq_Data/Samples/annotation_Groups_mCherry.csv
```

### 2. Add a New Line

Each row defines a new group using the following format:

| Group Name | Method           | Target Groups or Substrings         |
|------------|------------------|-------------------------------------|
| Non_target | include          | `['Liver', 'Muscle', 'Testis']`     |

### Supported Methods

- `include`: Combine listed groups.
- `exclude`: Include everything *except* the listed groups.
- `contains_include`: Include groups that **contain** a substring.
- `contains_exclude`: Include groups that **do not contain** a substring.

**Important**: The last column must be a Python list of strings. If there's more than one element, enclose it in square brackets to prevent column misalignment.

### 3. Re-run the Grouping Step

Once you've added your new group:

```bash
./Python_scripts/S6.py
```

This will update the final analysis tables with the new group included.

---

## How to Open the UI with Streamlit

An interactive web UI is available for exploring runs and results.

### 1. Activate the Conda Environment

The environment is already installed on the Monster machine, but you can recreate it using the YAML file in the `env/` folder:

```bash
conda env create -f env/report_env.yml
conda activate report
```

### 2. Run the Streamlit App

```bash
streamlit run Python_scripts/app.py
```

Then open your browser to view and interact with the analysis results.

---

## 🛠 Troubleshooting & Tips

- Always double-check your folder and file names — exact naming is required!
- When editing CSV files, ensure proper formatting (especially for list-style entries).
- If a step fails, look at the corresponding log file in the `logs/` directory for helpful error messages.
- If you have any other issues or questions, please reach out to me on Teams or via email.

---