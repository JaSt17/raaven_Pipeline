# Detailed Explanation of the Config Script

This config script holds **six dictionaries**, one for every step of the pipeline. This document provides detailed information on each parameter and how to set them.

---

## Global Parameters

There are two global parameters that **must** be set:

1. **data_dir**  
   - The directory where the entire pipeline is run. This is where input data is stored and where all results/logs will be saved.
2. **log_dir**  
   - By default: `data_dir + "/logs/"`
   - This is the folder where the logs from each step of the pipeline will be stored.
   - **Note**: Every dictionary for each step contains the path to `log_dir` to ensure the logfile is saved in the correct location.

---

## Config_S1

1. **input_file**  
   Path to the fasta file of all the reference proteins that were used to build the library.

2. **wSet**  
   - Path to the file containing the percentages for Homo sapiens codon usage values.
   - The file in `Example/input/` can be used as a default.
   - If replaced, the file must maintain the same format as the default, including the column starting with "hsa".

3. **structure_dict**  
   A dictionary containing the different structures of fragments that should be created from the reference file. Each structure has a name (used as the dictionary key) and the value is another dictionary with the following keys:
   - **length**: Length of the fragment in number of amino acids.
   - **freq**: Step size of the sliding window method. A value of `1` shifts the window by 1 amino acid (3 bases) at a time.
   - **overhangs**: A tuple of the overhangs that will be used for the fragments.

4. **output_csv**  
   The path where the Look-Up Table (LUT) will be saved.

5. **output_name**  
   The path where a list of all created fragments will be saved.

---

## Config_S2

1. **in_name_barcode**  
   Path to the fastq file of the library barcode reads.

2. **in_name_fragment**  
   Path to the fastq file of the library fragment reads.

3. **out_dir**  
   Path to the directory where we will save the found barcodes and fragments.

4. **out_name**  
   Name that will be appended to the found fragments and barcode files.

5. **bbduk2_args**  
   Two lists of arguments for the `bbduk2` runs. Below is a short description of each parameter:

   - **k**: Specifies the k-mer length used for matching.
   - **hammingdistance**: Number of allowed mismatches when comparing sequences.
   - **overwrite**: Enables overwriting of output files if they already exist.
   - **findbestmatch**: Searches for the best match among multiple possibilities.
   - **rcomp**: Enables reverse complement matching.
   - **qhdist**: Number of mismatches in quality-based Hamming distance for matching.
   - **minavgquality**: Minimum average quality requirement for sequences.
   - **maxns**: Number of `N` bases allowed in the sequences.
   - **minlength**: Minimum length of the sequences.
   - **maxlength**: Maximum length of the sequences.
   - **ordered**: Maintains the order of sequences in the output.
   - **lliteral**: Specifies the left adapter sequence to match.
   - **rliteral**: Specifies the right adapter sequence to match.

   > **Note**: The literals contain the overhangs from **S1**.

---

## Config_S3

1. **in_name_LUT, barcode_file, fragment_file**  
   - All set according to the names given in **Config_S1** and **Config_S2**.

2. **threshold**  
   - Minimal ratio of the most frequent barcode to all found barcodes for chimeric barcode detection.
   - If the ratio of a barcode matched to one fragment is higher than the threshold, it is classified as a valid barcode; if lower, it is classified as chimeric.

3. **chunk_size**  
   - Number of lines processed at a time from the barcode and fragment files.
   - Adjust based on the computing power of the system. Smaller values increase processing time but reduce memory requirements, and vice versa.

4. **out_name**  
   - Path to the output file that will be created at the end of **S3**.

---

## Config_S4

1. **input_table, in_name_LUT**  
   - Set according to the names given in **Config_S1** and **Config_S3**.

2. **sample_inputs**  
   - Path to the list containing the names of all tissue samples to be analyzed for barcodes.
   - The list should have two columns:  
     1. The name of the file.  
     2. The tissue or group from which we have the sequencing data.

3. **sample_directory**  
   - Path to the directory that contains the sequencing data for the tissues/samples listed in `sample_inputs`.

4. **log_file_path**  
   - Path to the file where the details of the found barcodes will be saved.

5. **output_dir**  
   - Path to the directory where all found barcode files will be saved.

6. **bbduk2_args**  
   - A list of arguments for the `bbduk2` runs to extract the barcodes.
   - Below are parameters not yet explained in **Config_S2**:
     - **trd**: Trims adapters from both the read and its reverse complement.
     - **skipr2**: Skips processing of Read 2 in paired-end data.

---

## Config_S5

1. **input_table, in_name_LUT**  
   - Set according to the names given in **Config_S1** and **Config_S3**.

2. **output_table**  
   - Path to the output file of **S5**.

---

## Config_S6

1. **original_seq_file, input_dir, sample_inputs, library_fragments**  
   - Set according to the names given in **Config_S1**, **Config_S4**, and **Config_S5**.

2. **library_name**  
   - Group name that will be assigned to all the fragments found in the Library sequencing data in the final output.

3. **subsets**  
   - A dictionary containing information about different subsets to be created in the final output.
   - The key is the name of the subset group, and the value is a list of the fragments that should be included or excluded.

4. **backbone_seq**  
   - The backbone sequence added to the fragments, which will be cut away in this last step.

5. **trim_dict**  
   - A dictionary holding the positions for trimming for every structure created in **S1**.
   - These trimming positions depend on the overhangs that were given to the sequences.

6. **output_table**  
   - Path to the final output file of the pipeline.

> **Note**: The `backbone_seq` and `trim_dict` do not have to be set. They are only necessary depending on how the fragments were cut during the `bbduk2` step in **S2**.

---
