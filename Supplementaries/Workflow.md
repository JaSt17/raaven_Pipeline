Python Pipeline for Extracting Barcodes & Fragments and Creating a Quantified Output Table
========================================================================================

Overview and Key Clarifications
-------------------------------
1. Literal Sequences:
   - Chosen based on the adapters used during library creation.
   - Serve as the anchors (“left” and “right” flanks) for identifying and extracting barcodes/fragments in raw sequencing reads.

2. Sequencing Setup:
   - Paired-End Illumina Reads (R1 and R2) without overlap.
   - Barcode is found on Read 1 (R1), while the fragment is on Read 2 (R2).

3. Software Environment:
   - bbduk2 version is specified in the README.
   - Other tools (e.g., starcode, seqkit) and Python packages are listed in env/pipeline.yml.

4. Distance & Threshold Configurations:
   - Hamming distances for matching literals (±2) and barcodes (±1) are configurable in the pipeline’s config files.
   - Ratio threshold for allowing chimeric barcodes is configurable (currently set to 100% – meaning no chimeric allowance).


Step 1: Fragment Generation and Translation
---------------------------------------------
Purpose:
   - Generate all possible DNA inserts (fragments) from reference protein sequences using a sliding window approach.
   - Record relevant metadata in a lookup table (LUT).

Processes:
   1. Translation & Back-Translation:
      - Translation: Convert original DNA sequences to amino acids.
      - Back-Translation: Convert amino acids back to DNA using a Homo sapiens–optimized codon table.
      - Output: A reference DNA sequence file (with human-optimized codons) used later for validating fragment identity.

   2. Fragmentation:
      - Use a sliding window approach with a defined fragment insertion size to generate all possible fragments.
      - Accommodate multiple structures if fragment overhangs vary.

   3. Lookup Table (LUT) Creation:
      - Store each fragment in a structured dataframe with metadata:
         * Protein name of origin.
         * Fragment structure type (e.g., different overhang designs).
         * Position within the protein.
         * Peptide sequence (amino acid translation).

Output of Step 1:
   - A comprehensive LUT detailing all fragments, their metadata, and references for downstream matching.


Step 2: Barcode and Fragment Extraction
-----------------------------------------
Purpose:
   - Isolate barcodes and fragments from raw Illumina FASTQ files using known adapter-based “literal” sequences.

Processes:
   1. Defining the Literals:
      - Left/Right Flanks: Based on the adapters used during library creation.
      - Configurable Hamming Distance: Set to 1 to allow minor sequencing errors.

   2. Fragment Extraction (R2):
      - Target: Fragments typically appear on Read 2.
      - Using bbduk2:
         * Extract all potential fragment reads using the literal.
         * Retain only reads with an exact fragment length.
      - Matching to Reference Sequences (if a reference exists):
         * Using VSearch, compare extracted fragments to the reference DNA (from Step 1) with 100% identity.
         * Discard fragments that fail to match exactly.
         * If no reference is provided, use all extracted fragments.
      - Extract matching fragments and corresponding barcode reads:
         * Using Seqkit, extract matching reads from both files based on headers identified by VSearch.

   3. Barcode Extraction (R1):
      - Target: Barcodes typically appear on Read 1.
      - Using bbduk2:
         * Extract all correct barcodes using the literal.
         * Only consider reads that have a matching corresponding fragment.

   4. Pairing Barcodes with Fragments:
      - Use seqkit pair to match each barcode (R1) with its corresponding fragment (R2) by read ID.
      - Store the paired data in new FASTQ files.

Output of Step 2:
   - Barcode–fragment pairs that will feed into Step 3 for clustering and classification.


Step 3: Barcode Reduction and Classification
-----------------------------------------------
Purpose:
   - Create a comprehensive Barcode–Fragment table, cluster similar barcodes to correct sequencing errors, and exclude ambiguous entries.

Processes:
   1. Extracting Unique Fragments & Barcodes:
      - Gather all unique barcodes and fragments as a reference for downstream steps.
      - Visualize library diversity.

   2. Building the Full Table:
      - Merge barcodes and fragments (possibly in chunks due to file size).
      - Each row links a fragment to a barcode and an entry in the LUT.

   3. Barcode Clustering (using Starcode):
      - Cluster barcodes within a Hamming distance of 1 (configurable).
      - Each cluster is represented by a consensus barcode, reducing duplicates caused by sequencing errors.
      - Replace original barcodes in the table with their cluster consensus versions.

   4. Filtering Barcodes:
      - Single-Read Barcodes:
         * Observed only once → Mark as single.
         * May be excluded from further analysis (configurable).
      - Multi-Read Barcodes:
         * May map to exactly one fragment (clean) or multiple fragments (chimeric).
         * For chimeric barcodes, calculate the ratio mCount/tCount for the most frequent fragment.
         * Current Threshold = 80%; barcodes that do not exclusively map to one fragment are deemed chimeric and excluded.

   5. Generating the Final Table:
      - Save valid barcode–fragment pairs in library_barcodes.csv.
      - Save single or chimeric barcodes separately (if needed).
      - Write a Barcode database file containing all valid barcodes from the Plasmid library, which is used to extract valid barcodes from tissue samples in the next step.

Output of Step 3:
   - A cleaned and consolidated table linking high-confidence barcodes to fragments.


Step 4: Sample Barcode Processing
----------------------------------
Purpose:
   - Extract and process barcodes from RNA samples using a similar approach to Step 2, applied to experimental data.

Processes:
   1. Loading Data:
      - Library Fragments (from Steps 1–3)
      - Lookup Table (LUT) with detailed metadata.
      - Sample Metadata: Paths, group names, etc.

   2. Barcode Extraction (R1):
      - Using bbduk2:
         * Extract all correct barcodes using the literal.

   3. Aligning Barcodes:
      - Align extracted barcodes to the Barcode Database from S3.
      - Using VSearch:
         * Compare extracted barcodes to the Barcode Database (from Step 3) with 95% identity.
         * Discard any barcode that fails to match.

   4. Matching Barcodes with Fragments:
      - Match barcodes to library fragments from Steps 1–3.
      - Append peptide sequences and relevant metadata.
      - Sort results by read frequency (RNA counts) and save per sample.

   5. Summary Statistics:
      - Log the number of reads, barcodes, consensus clusters, matched fragments, etc., for each sample.

Output of Step 4:
   - Per-sample CSV files listing matched fragments, read frequencies (RNA counts), and starcode cluster barcodes.


Step 5: Data Merging (Library Fragments + Reference Positions)
--------------------------------------------------------------
Purpose:
   - Merge library fragment data (from Steps 1–4) with additional positional information to produce an annotated table.

Processes:
   1. Loading Data:
      - Library Fragments: With barcode counts, modes, etc.
      - Fragments Position File: Contains start/end positions, widths, or other structural info.
      - Note: If the positional file is missing (e.g., in NNK libraries), only library fragment counts are stored.

   2. Merging & Annotation:
      - Merge based on LUT number (LUTnr) and/or peptide.
      - Remove unneeded columns, rename tCount to RNAcount for clarity, and reorder columns.

   3. Saving Annotated Data:
      - The output CSV contains:
         * Fragment Info: origin_seq, mode, structure, LUTnr, peptide.
         * Positional Details: start, end, width, sequence.
         * Counts: mCount, RNAcount.

Output of Step 5:
   - An annotated library fragments CSV, enriched with positional data.


Step 6: Final Dataset Processing and Normalization
---------------------------------------------------
Purpose:
   - Combine data from multiple samples, normalize read counts to account for sequencing depth, and optionally trim overhangs.

Processes:
   1. Loading & Combining Data:
      - Load CSVs from each sample directory.
      - Rename groups for consistency (using sample metadata).
      - Integrate library fragment data for completeness.
      - Optionally add sequence length info from a FASTA file.

   2. Subsetting Data:
      - Create subsets based on user-defined conditions (e.g., group membership).
      - Tag subsets and merge them back into the main dataset for comparative analysis.

   3. Normalizing Read Counts:
      - Perform group-wise normalization to enable fair comparisons across samples with different sequencing depths.

   4. Aggregating Identical Fragments:
      - Combine counts for identical fragments across groups.
      - Aggregate total reads (tCount, mCount, RNAcount, normalized counts).
      - Track unique barcodes (BC) and LUT numbers (LUTnr).
      - Optionally calculate a “barcode-adjusted count ratio” for fragments with multiple barcodes.

   5. Saving the Final Output:
      - Save the merged dataset as a CSV file, ready for statistical analysis or machine learning workflows.

Output of Step 6:
   - Final annotated dataset of all found barcodes with their corresponding fragments from all samples and the library.


Conclusion
----------
This pipeline guides you from reference protein sequences to high-confidence barcode–fragment associations. It processes RNA sample data and ultimately produces a normalized, annotated dataset. By customizing parameters such as Hamming distances and ratio thresholds in the config files, the workflow can be adapted to various sequencing platforms and quality requirements.
