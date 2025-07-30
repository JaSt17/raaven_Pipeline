# 📘 Pipeline Configuration README

This document describes the configuration parameters used in the pipeline. The configuration file is modular and divided into several steps (`S1` to `S6`), each corresponding to a stage in the pipeline.

---

## 🔧 Global Configuration Parameters

These are set at the top of the config file and are used across multiple pipeline steps.

```python
data_dir = "Projects/Kingfischer/Seq_Data"
save_dir = "Projects/Kingfischer/Libraries"
annotation_file = data_dir + "/Samples/annotation_p007.csv"
Group_annotation_file = data_dir + "/Samples/annotation_Groups.csv"
```

- `data_dir`: Path to raw sequencing data. There should be a subdirectory for each plasmid sequencing run, which corresponds to the `run` name given below. In that subdirectory, you should have the `<library_name>_R1` and `_R2` files for the fragment and barcode sequences. For full plasmid sequencing, you typically only need one `fastq.gz` file.
- `save_dir`: Path where all outputs (e.g., logs, plots, CSVs) will be stored. For each run, a new directory will be created in this path named after the `library_name`.
- `annotation_file`: Path to the annotation CSV file containing metadata for the mRNA samples. This file should be located in the `Samples` subdirectory of `data_dir`.
- `Group_annotation_file`: Path to the annotation CSV file containing metadata on how we want to create subgroups. This file should also be located in the `Samples` subdirectory of `data_dir`.

```python
library_name = "p007"
run = "Plasmid_run_1"
```

- `library_name`: Identifier for the library (e.g., `"p007"` from the Kingfischer project). This must match the file prefix of the `fastq.gz` files in the `Seq_Data` directory.
- `run`: Describes the sequencing run (e.g., `"Plasmid_run_1"`). This must match the subdirectory name under `Seq_Data`.

```python
bc_len = 27
frag_len = 27
linker_length = 3
```

- `bc_len`: Expected barcode length in bases.
- `frag_len`: Expected fragment length in bases plus the linkers. (21 pb fragment + 3 bp L-linker + 3 bp R-linker)
- `linker_length`: Number of bases flanking the fragment on both sides. Usually 3, as there's one alanine to the left and right of the insert.

```python
num_possible_frag = 117126
```

- This number can either be set after running `S1` (by counting how many reference sequences were generated), or manually if you already know the intended number of fragments.

```python
barcode_left_literal = "TATCGCAAGA"
barcode_right_literal = "ATAACTTCGT"
fragment_left_literal = "GAGAGGCAACG"
fragment_right_literal = "AGACAAGCAG"
```

- These are the literal sequences expected to flank the barcode and fragment regions. Used by `bbduk2` for accurate extraction.

```python
single_read_barcodes = True
chimeric_barcodes = False
starcode_reduction = False
threshold = 1
```

- `single_read_barcodes`: Allows fragment-barcode pairs that appear only once.
- `chimeric_barcodes`: Enables inclusion of chimeric barcodes.
- `starcode_reduction`: Clusters similar barcodes using a 1-base Hamming distance. (This most of the time increase the number of chimeric reads)
- `threshold`: Threshold for accepting chimeric reads with a dominant barcode-fragment pairing. For Illumina, keep it disabled. For PacBio, you could use `0.8` to accept barcodes where the most frequent fragment is 8× more common than the next.

---

## 🧬 Step S1: Create Fragment Lookup Table (LUT)

```python
config_S1 = {
    "input_file": ...,
    "structure_dict": {
        "7aa": {"length": 7, "freq": 1}},
    "LibID": ...,
    "output_csv": ...,
    "log_dir": ...
}
```

- `input_file`: FASTA file with reference protein sequences. This should be located in an `input/` folder within the save directory.
- `structure_dict`: Defines peptide fragmenting logic. `length` is amino acids per fragment, `freq` is the sliding window step (e.g., `1` = every amino acid is a start).
- All other parameters are set automatically from global variables and typically don’t need to be changed.

---

## 🧬 Step S2: Extract Barcode and Fragment Sequences

```python
config_S2 = {
    "in_name_barcode": ...,
    "in_name_fragment": ...,
    "draw_sequence_logos": True,
    ...
    "bbduk2_args_BC": [...],
    "bbduk2_args_Frag": [...],
}
```

- `draw_sequence_logos`: If `True`, generates sequence logo plots from input reads. Useful for visualizing structure in Illumina reads. Not suitable for PacBio reads due to length (>5000 bp).
- `bbduk2_args_BC` / `bbduk2_args_Frag`: Parameters used by `bbduk2` to extract barcodes and fragments. If you modify literal sequence lengths, update the `k` value to match.
- All other parameters are inherited from the global configuration and usually do not require changes.

---

## 🧬 Step S3: Create Library Barcode Table

```python
config_S3 = {
    "in_name_LUT": ...,
    "barcode_file": ...,
    "fragment_file": ...,
    ...
    "chunk_size": 20000000,
    "out_name": ...,
}
```

- `chunk_size`: Defines how many reads are processed at once. For example, on Lunarc you can use `20000000`. If you're running on a machine with less memory, consider reducing this to avoid OOM (out-of-memory) errors.
- All other parameters are automatically inherited and do not need to be modified.

---

## 🧬 Step S4: Match Barcodes in Experimental Samples

```python
config_S4 = {
    "input_table": ...,
    "sample_inputs": ...,
    "in_name_LUT": None,
    "sample_directory": ...,
    "log_file_path": ...,
    ...
    "bbduk2_args": [...],
}
```

- `in_name_LUT`: If set, only barcodes present in the library reference are matched. Set to `None` to allow detection of new/unknown barcodes.
- `sample_inputs`: here you need a path to a annotation.csv file which holds the file names of all files that you want to extract barcodes from as well as their Tissue type.
- `sample_directory`: This should be the path tho the RNA sequenceing data.
- `bbduk2_args`: Extraction parameters for matching barcodes in samples.
- All other parameters are inherited from the global configuration and typically don't need to be changed.

---

## 🧬 Step S5: Map Found Barcodes to Fragment Information

```python
config_S5 = {
    "input_table": ...,
    "output_table": ...,
    ...
}
```

- All parameters are auto-generated from previous steps or global settings. No need for manual changes here.

---

## 🧬 Step S6: Generate Summary and Subset Data

```python
config_S6 = {
    "original_seq_file": ...,
    "LUT_file": ...,
    "input_dir": ...,
    "library_fragments": ...,
    "annotation_groups": ...,
    ...
    "output_table": ...,
```

- All parameters are auto-generated from previous steps or global settings. No need for manual changes here.

## 🔄 Configuration Access Function

```python
def get_config(step: str) -> dict:
    ...
```

- Simple helper to fetch the config dictionary for a given step like `"S1"` or `"S4"`.

---
