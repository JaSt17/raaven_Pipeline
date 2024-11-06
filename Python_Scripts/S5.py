#!/usr/bin/env python3
"""
Barcoded extraction and reduction from all given Samples.

Author: Jaro Steindorff

This script extracts barcodes from given samples, reduces them using the Starcode algorithm, and identifies
corresponding fragments. It processes a csv file with file_path and basename and saves the results for each base name group.

Workflow:
    - Load necessary data from previous steps (LUT.csv, MatchedFragments.csv, fragment_pos.csv)
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
import re
import logging
from datetime import datetime
import pandas as pd
from multiprocessing import Pool, cpu_count
import gzip
from Bio import SeqIO
from functools import partial
# local import
from config import get_config

# Initialize logging with custom format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%H:%M:%S',  # Only show hour, minute, and second
    filemode='w',  # Overwrite log file
    filename='Python_Scripts/Logs/S5.log'  # Log file name
    )
logger = logging.getLogger(__name__)

def run_command(command, description: str, shell=False, verbose=False):
    """
    Runs a subprocess command and returns stdout, stderr.

    :param command: Command list to execute or string if shell=True
    :param description: Description of the command for logging
    :param shell: Whether to execute the command through the shell
    """
    if verbose:
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
        
        
# Function to index fragments that are found multiple times
def match_range(idx_frag: int, foundFrags: pd.DataFrame, all_fragments_ranges : pd.DataFrame)->list:
    """ Helper function to find the ranges of a fragment in all_fragments_ranges based on its index in foundFrags.
    
    :param idx_frag: Index of the fragment in foundFrags
    :param foundFrags: DataFrame with the found fragments
    :param all_fragments_ranges: DataFrame with the ranges of all fragments
    
    :return: List of tuples with the index of the fragment in all_fragments_ranges and the index of the fragment in foundFrags
    """
    fragment_seq = foundFrags.iloc[idx_frag]['fragment']
    match_ranges = all_fragments_ranges[all_fragments_ranges['Sequence'] == fragment_seq].index.tolist()
        
    return [(match_range, idx_frag) for match_range in match_ranges]


def analyze_tissue(file_path:str, base_name:str, data_dir:str, out_dir:str, library_fragments: pd.DataFrame,
                    all_fragments_ranges: pd.DataFrame, lut_dna: pd.DataFrame, threads:int, bbduk2_args: list) -> dict:
    """
    Analyze a single tissue sample based on its index in the load list.

    :param file_path: path to the fastq file
    :param base_name: base name of the file
    :param data_dir: Directory containing the fastq files
    :param library_fragments: DataFrame loaded from multipleContfragmentsComplete.pkl
    :param all_fragments_ranges: DataFrame loaded from alignedLibraries.pkl
    :param lut_dna: DataFrame loaded from LUTdna.pkl
    :param threads: Number of threads to use
    :return: Dictionary containing log information for the sample
    """
    log_entry = {}
    # Add the name to the log entry and set the output name
    log_entry['Name'] = base_name
    path = os.path.join(data_dir, file_path)

    # Count the reads in the FASTQ file
    try:
        with gzip.open(path, "rt") as handle:
            read_count = sum(1 for _ in SeqIO.parse(handle, "fastq"))
    except Exception as e:
        logger.info(f"Error reading FASTQ file {path}: {e}")
        return None
    log_entry['Reads'] = read_count
    name_out = base_name
    
    logger.info(f"Processing {file_path} with {read_count} reads")

    # Extraction of barcodes
    # ============================
    # Create a temporary file for the output
    out_name_BC = tempfile.NamedTemporaryFile(prefix="BC_", suffix=".fastq.gz", delete=False).name

    # Run the bbduk2 command to extract barcodes from the reads
    bbduk2_command = ["bbmap/bbduk2.sh"] + bbduk2_args + [
        f"threads={threads}",
        f"in={path}",
        f"out={out_name_BC}",
    ]
    stdout, stderr = run_command(bbduk2_command, f"bbduk2 barcode extraction for {file_path}")

    # Save the Barcode extraction result in the log entry
    match = re.search(r"Result:\s*(\d+)", stderr)
    if match:
        log_entry['BCs'] = int(match.group(1))
    else:
        log_entry['BCs'] = 0

    # Read the barcodes
    try:
        with gzip.open(out_name_BC, "rt") as handle:
            reads_BC = list(SeqIO.parse(handle, "fastq"))
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
    _, stderr = run_command(starcode_command, f"starcode barcode reduction for {file_path}", shell=True)

    # Create a table with the reduced barcodes
    try:
        table_BC_sc = pd.read_csv(out_name_BC_star, sep="\t", header=None, names=['scBC', '_', 'BC'], dtype=str)
    except Exception as e:
        logger.error(f"Error reading Starcode output file {out_name_BC_star}: {e}")
        return None
    # Split the barcodes by comma
    table_BC_sc['BC'] = table_BC_sc['BC'].str.split(',')
    # Explode the table to have one row per barcode
    table_BC_sc = table_BC_sc.explode('BC')

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
    
    # Matching barcodes with fragments and adding RNAcount
    # ============================
    # Reorganize the barcode_table to have BC and RNAcount
    BCcount = barcode_table['BC'].value_counts().reset_index()
    BCcount.columns = ['BC', 'RNAcount']
    # Extract only BC that are in BCcount
    foundFrags = library_fragments.merge(BCcount, on='BC', how='inner')
    # Merge with lut_dna on 'LUTnr'
    foundFrags = foundFrags.merge(lut_dna, on='LUTnr', how='inner')
    # Rename 'Sequence' to 'fragment'
    foundFrags.rename(columns={'Sequence': 'fragment'}, inplace=True)
    # Delete the 'Names' and 'i.Structure' columns if they exist
    for col in ['Names', 'i.Structure']:
        if col in foundFrags.columns:
            foundFrags.drop(columns=[col], inplace=True)

    # Create a list of matches
    # create a pool of workers
    pool = Pool(threads)
    # Create a partial function to pass the necessary arguments
    match_range_partial = partial(match_range, foundFrags=foundFrags, all_fragments_ranges=all_fragments_ranges)
    # Match the ranges of the fragments by using all available threads
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
        output_filename = os.path.join(out_dir, f"found.{name_out}.csv")
        found_fragments_ranges.to_csv(output_filename, index=False)
    else:
        logger.info(f"No matches found for sample {name_out}")
    
    logger.info(f"Finished processing {file_path} found {log_entry['scBCs']} barcodes")    
    
    return log_entry

def main():
    start_time = datetime.now()
    threads = cpu_count()
    
    # load config
    config = get_config("S5")

    # Try to load the necessary data
    try:
        library_fragments = pd.read_csv(config["input_table"], dtype={7: str})
        all_fragments_ranges = pd.read_csv(config["fragments_pos"])
        lut_dna = pd.read_csv(config["in_name_LUT"])
    except Exception as e:
        logger.error(f"Error loading data: {e}")
        sys.exit(1)
    # Try to load the sample inputs
    try:
        load_list = pd.read_csv(config["sample_inputs"], header=None)
    except Exception as e:
        logger.error(f"Error loading sample inputs: {e}")
        sys.exit(1)
    # check if the data directory exists
    data_dir = config["sample_directory"]
    output_dir = config["output_dir"]
    # Create the output directory if it does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(data_dir):
        logger.error(f"Data directory {data_dir} does not exist")
        sys.exit(1)
    log_table = []
    # get the settings for the barcode extraction
    bbduk2_args_BC = config["bbduk2_args"]

    # Analyze each tissue sample
    for row in load_list.iterrows():
        # Extract the file name from the first column
        file_path = row[1][0]
        # Extract the base name from the second column
        base_name = row[1][1]
        log_entry = analyze_tissue(file_path, base_name, data_dir, output_dir, library_fragments, all_fragments_ranges, lut_dna, threads, bbduk2_args_BC)
        if log_entry:
            log_table.append(log_entry)

    # Create a DataFrame from the log table
    log_df = pd.DataFrame(log_table)
    # Save the log table
    log_df.to_csv(config["log_file_path"], index=False)

    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
