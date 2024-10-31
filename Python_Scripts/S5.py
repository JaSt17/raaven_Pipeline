#!/usr/bin/env python3
"""
Barcoded extraction and reduction from RNA samples

Author: Jaro Steindorff

This script extracts barcodes from RNA samples, reduces them using the Starcode algorithm, and identifies
corresponding fragments. It processes a list of RNA samples and saves the results for each sample.

Workflow:
    - Load necessary data from previous steps
    - For each RNA sample:
        - Extract barcodes using bbduk2.sh
        - Reduce barcodes using Starcode
        - Match reduced barcodes with fragments
        - Save found fragments for the sample
    - Save a log table with summary statistics
"""

import os
import sys
import subprocess
import tempfile
import shutil
import re
import logging
from datetime import datetime
import pandas as pd
from multiprocessing import Pool, cpu_count
from Bio import SeqIO
from functools import partial
# local import
from config import get_config

# Initialize logging with custom format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%H:%M:%S'  # Only show hour, minute, and second
)
logger = logging.getLogger(__name__)

def run_command(command, description: str, shell=False):
    """
    Runs a subprocess command and returns stdout, stderr.

    :param command: Command list to execute or string if shell=True
    :param description: Description of the command for logging
    :param shell: Whether to execute the command through the shell
    """
    logger.info(f"Running {description}")
    try:
        if shell:
            process = subprocess.run(command, shell=True, capture_output=True, text=True)
        else:
            process = subprocess.run(command, capture_output=True, text=True)
        stdout = process.stdout
        stderr = process.stderr
        if process.returncode != 0:
            logger.error(f"Error running {description} with code {process.returncode}")
            logger.error(stderr)
            sys.exit(1)
        return stdout, stderr
    except Exception as e:
        logger.error(f"Exception running {description}: {e}")
        sys.exit(1)

def analyze_tissue(index_nr, load_list, data_dir, output_table, all_fragments_ranges, lut_dna, threads):
    """
    Analyze a single tissue sample based on its index in the load list.

    :param index_nr: Index number of the sample in the load list
    :param load_list: DataFrame containing the load list
    :param data_dir: Directory containing the fastq files
    :param output_table: DataFrame loaded from multipleContfragmentsComplete.pkl
    :param all_fragments_ranges: DataFrame loaded from alignedLibraries.pkl
    :param lut_dna: DataFrame loaded from LUTdna.pkl
    :param threads: Number of threads to use
    :return: Dictionary containing log information for the sample
    """
    log_entry = {}
    # Extract the base name for the file corresponding to the provided index
    base_name = load_list.loc[index_nr, 'Name']

    # Split the base name into components using '/' as a delimiter
    name_components = [comp for comp in base_name.split('/') if comp]

    # Determine the input file paths based on the number of name components
    if len(name_components) == 2:
        # If the name has two components, construct the path to the files within the subdirectory
        path = os.path.join(data_dir, name_components[0])
        pattern = name_components[1]
    else:
        # If the name has a single component, search for files directly within data_dir
        path = data_dir
        pattern = name_components[0]

    # Get list of files matching the pattern
    try:
        in_files = [os.path.join(path, f) for f in os.listdir(path) if f.startswith(pattern)]
    except FileNotFoundError:
        logger.error(f"Path not found: {path}")
        return None

    # Identify the files corresponding to P5 (e.g., '_1' in the name) because the barcodes are in the P5 reads
    in_name_P5 = [f for f in in_files if '_1' in f]
    if not in_name_P5:
        logger.error(f"No P5 files found for index {index_nr}")
        return None
    in_name_P5 = in_name_P5[0]

    # Add the name to the log entry and set the output name
    log_entry['Name'] = base_name

    # Count the number of reads in in_name_P5
    try:
        read_count = sum(1 for _ in SeqIO.parse(in_name_P5, "fastq"))
    except Exception as e:
        logger.error(f"Error reading FASTQ file {in_name_P5}: {e}")
        return None
    log_entry['Reads'] = read_count
    name_out = base_name

    # Extraction of barcodes
    # ============================
    # Create a temporary file for the output
    out_name_BC = tempfile.NamedTemporaryFile(prefix="BC_", suffix=".fastq.gz", delete=False).name

    # Run the bbduk2 command to extract barcodes from the reads
    bbduk2_command = [
        "bbmap/bbduk2.sh",
        "overwrite=true",
        "k=12",
        "mink=12",
        "hammingdistance=2",
        "findbestmatch=t",
        "trd=t",
        "rcomp=f",
        "skipr2=t",
        "findbestmatch=f",
        "qhdist=0",
        "minavgquality=0",
        "ordered=t",
        "maxns=0",
        "minlength=18",
        "maxlength=22",
        f"threads={threads}",
        f"in={in_name_P5}",
        f"out={out_name_BC}",
        "lliteral=GGCCTAGCGGCCGCTTTACTT",
        "rliteral=ATAACTTCGTATA"
    ]
    stdout, stderr = run_command(bbduk2_command, "bbduk2 barcode extraction")

    # Save the Barcode extraction result in the log entry
    match = re.search(r"Result:\s*(\d+)", stderr)
    if match:
        log_entry['BCs'] = int(match.group(1))
    else:
        log_entry['BCs'] = 0

    # Read the barcodes
    try:
        reads_BC = list(SeqIO.parse(out_name_BC, "fastq"))
    except Exception as e:
        logger.error(f"Error reading FASTQ file {out_name_BC}: {e}")
        return None
    # Create a table with the barcodes and the ID
    ids = [record.id for record in reads_BC]
    bcs = [str(record.seq) for record in reads_BC]
    barcode_table = pd.DataFrame({'ID': ids, 'BC': bcs})

    # Starcode based barcode reduction
    # ============================
    # Create a temporary file for the output
    out_name_BC_star = tempfile.NamedTemporaryFile(prefix="BCsc_", suffix=".txt", delete=False).name
    # Run the starcode command to reduce the barcodes
    starcode_command = f"gunzip -c {out_name_BC} | starcode -t {threads - 1} --print-clusters -d 1 -r5 -q -o {out_name_BC_star}"
    _, stderr = run_command(starcode_command, "starcode barcode reduction", shell=True)

    # Create a table with the reduced barcodes
    try:
        table_BC_sc = pd.read_csv(out_name_BC_star, sep="\t", header=None, index_col=0, names=['V1', 'V2', 'V3'], dtype=str)
    except Exception as e:
        logger.error(f"Error reading Starcode output file {out_name_BC_star}: {e}")
        return None
    table_BC_sc.reset_index(inplace=True)
    table_BC_sc.rename(columns={'index': 'rn'}, inplace=True)
    # Split the third column by comma
    table_BC_sc['V3'] = table_BC_sc['V3'].str.split(',')
    table_BC_sc = table_BC_sc.explode('V3')
    # Set the column names to BC and scBC
    table_BC_sc.rename(columns={'rn': 'scBC', 'V3': 'BC'}, inplace=True)

    # Replacing barcodes with Starcode reduced versions
    # ============================
    # Merge barcode_table with table_BC_sc on 'BC'
    barcode_table = barcode_table.merge(table_BC_sc[['BC', 'scBC']], on='BC', how='inner')
    # Rename columns
    barcode_table.rename(columns={'BC': 'oldBC', 'scBC': 'BC'}, inplace=True)
    # Save the number of unique old barcodes and remaining barcodes in the log entry
    all_BCs = barcode_table['oldBC'].nunique()
    sc_BCs = barcode_table['BC'].nunique()
    SC_dropped_BC = all_BCs - sc_BCs
    log_entry['allBCs'] = all_BCs
    log_entry['scBCs'] = sc_BCs
    log_entry['SCdroppedBC'] = SC_dropped_BC
    # Delete the old barcodes column
    barcode_table.drop(columns=['oldBC'], inplace=True)

    # Reorganize the barcode_table to have BC and RNAcount
    BCcount = barcode_table['BC'].value_counts().reset_index()
    BCcount.columns = ['BC', 'RNAcount']
    # Extract only BC that are in BCcount
    foundFrags = output_table.merge(BCcount, on='BC', how='inner')
    # Merge with lut_dna on 'LUTnr'
    foundFrags = foundFrags.merge(lut_dna, on='LUTnr', how='inner')
    # Rename 'Sequence' to 'fragment'
    foundFrags.rename(columns={'Sequence': 'fragment'}, inplace=True)
    # Delete the 'Names' and 'i.Structure' columns if they exist
    for col in ['Names', 'i.Structure']:
        if col in foundFrags.columns:
            foundFrags.drop(columns=[col], inplace=True)

    # Function to index fragments that are found multiple times
    def match_range(idx_frag, foundFrags, all_fragments_ranges):
        fragment_seq = foundFrags.iloc[idx_frag]['fragment']
        match_ranges = all_fragments_ranges[all_fragments_ranges['Sequence'] == fragment_seq].index.tolist()
        return [(match_range, idx_frag) for match_range in match_ranges]

    # Create a list of matches
    pool = Pool(threads)
    match_range_partial = partial(match_range, foundFrags=foundFrags, all_fragments_ranges=all_fragments_ranges)
    match_ranges_list = pool.map(match_range_partial, range(len(foundFrags)))
    pool.close()
    pool.join()
    # Flatten the list of matches
    match_ranges = [item for sublist in match_ranges_list for item in sublist]

    if match_ranges:
        match_ranges_df = pd.DataFrame(match_ranges, columns=['matchRanges', 'idxFrag'])
        # Select foundFrags based on idxFrag
        foundFrags = foundFrags.iloc[match_ranges_df['idxFrag']].reset_index(drop=True)
        # Remove unnecessary columns
        for col in ['Reads', 'fragment', 'Structure', 'LUTnr']:
            if col in foundFrags.columns:
                foundFrags.drop(columns=[col], inplace=True)
        # Select all_fragments_ranges based on matchRanges
        found_fragments_ranges = all_fragments_ranges.iloc[match_ranges_df['matchRanges']].reset_index(drop=True)
        # Add foundFrags columns to found_fragments_ranges
        found_fragments_ranges = pd.concat([found_fragments_ranges, foundFrags], axis=1)
        # Sort them by RNAcount in descending order
        found_fragments_ranges.sort_values(by='RNAcount', ascending=False, inplace=True)
        # Save the found fragments for the sample in the output folder
        output_filename = os.path.join("2_output", f"found.{name_out}.pkl")
        os.makedirs(os.path.dirname(output_filename), exist_ok=True)
        found_fragments_ranges.to_pickle(output_filename)
    else:
        logger.info(f"No matches found for sample {name_out}")
    return log_entry

def main():
    start_time = datetime.now()
    threads = cpu_count()
    
    # load config
    config = get_config(5)

    # Load the data
    try:
        output_table = pd.read_csv(config["input_table"])
        all_fragments_ranges = pd.read_csv(config["fragments_pos"])
        lut_dna = pd.read_csv(config["in_name_LUT"])
        load_list = pd.read_csv(config["sample_inputs"], sep="\t", header=None, names=["Name", "BaseName", "GroupName"])
    except Exception as e:
        logger.error(f"Error loading data: {e}")
        sys.exit(1)

    data_dir = config["sample_directory"]
    log_table = []

    # Analyze each tissue sample
    for index_nr in range(len(load_list)):
        log_entry = analyze_tissue(index_nr, load_list, data_dir, output_table, all_fragments_ranges, lut_dna, threads)
        if log_entry:
            log_table.append(log_entry)

    # Create a DataFrame from the log table
    log_df = pd.DataFrame(log_table)
    # Save the log table
    log_df.to_pickle("S5_log.table.pkl")

    # Cleanup the temp files
    temp_dir = tempfile.gettempdir()
    try:
        shutil.rmtree(temp_dir)
    except Exception as e:
        logger.error(f"Error cleaning up temp files: {e}")

    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
