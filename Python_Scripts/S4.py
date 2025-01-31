#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script extracts barcodes from given samples, reduces them using the Starcode algorithm, and matches them with
corresponding fragments. It processes a csv file with Sample and Group and saves the results in a log table and saves the found fragments in a csv file.

Workflow:
    - Load the library fragments and LUT data
    - For each RNA sample:
        - Extrating known barcodes from the sample using vsearch
        - Extract barcodes using bbduk2.sh
        - Reduce barcodes using Starcode
        - Match reduced barcodes with fragment information
        - Save found fragments for the sample
    - Save a log table with summary statistics
    
Input:
    - library_fragments: DataFrame with the library fragments
    - fragments_pos: DataFrame with the positions of the fragments in the original sequneces
    - lut_dna: DataFrame with the LUT data
    - sample_inputs: CSV file with the sample inputs (file_path, base_name)
    - sample_directory: Directory containing the fastq files
    - log_file_path: Path to save the log table
    - output_dir: Directory to save the found fragments
    - bbduk2_args: List of arguments for bbduk2.sh barcode extraction
    
"""

import os
import sys
import subprocess
import tempfile
import re
import shutil
import logging
from datetime import datetime
import pandas as pd
import multiprocessing
import gzip
import gc
from Bio import SeqIO
from itertools import islice
# local import
from config import get_config

# function to create a global logger
def create_logger(path: str, name: str) -> None:
    """
    Create a global logger with a custom format.
    
    Parameters:
        path (str): The path to the log file
        name (str): The name of the logger
        
    Returns:
        None
    """
    filename = path + name + ".log"
    # Initialize logging with custom format
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        datefmt='%H:%M:%S',  # Only show hour, minute, and second
        filemode='w',  # Overwrite log file
        filename=filename
        
    )
    global logger  # Declare the logger as global
    logger = logging.getLogger(name) # Create a logger
    

def run_command(command, description: str, shell=False, verbose=False):
    """
    Runs a subprocess command and returns stdout, stderr, and error status.
    
    Parameters:
        command (list): The command to run
        description (str): A description of the command
        
    Returns:
        tuple: The stdout and stderr of the command
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
        
        
def starcode_based_barcode_reduction(full_table: pd.DataFrame) -> pd.DataFrame:
    """
    Perform barcode reduction using Starcode clustering on the unique barcodes from the full table.

    Parameters:
        full_table (pd.DataFrame): DataFrame containing the full table with barcodes (column 'BC')
        
    Returns:
        pd.DataFrame: DataFrame containing the Starcode-reduced barcodes
    """
    
    # Write unique barcodes to a temporary file
    barcode_temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w+')
    barcode_temp_file.writelines(f"{bc}\n" for bc in full_table['BC'])
    barcode_temp_file.close()
    
    # Run Starcode clustering
    num_threads = multiprocessing.cpu_count()
    starcode_output = tempfile.NamedTemporaryFile(delete=False, suffix=".txt")
    starcode_cmd = f"starcode -t {num_threads-1} --print-clusters -d 1 -r5 -q -i {barcode_temp_file.name} -o {starcode_output.name}"
    subprocess.run(starcode_cmd, shell=True, check=True)

    # Read Starcode output
    starcode_df = pd.read_csv(starcode_output.name, sep='\t', header=None, names=['scBC', 'count', 'BC_list'])
    starcode_df['BC_list'] = starcode_df['BC_list'].str.split(',')

    # Explode the BC_list to create a mapping between original barcodes and Starcode-reduced barcodes
    starcode_exploded = starcode_df.explode('BC_list')
    starcode_exploded.rename(columns={'BC_list': 'BC'}, inplace=True)

    # Clean up temporary files
    for temp_file in [barcode_temp_file.name, starcode_output.name]:
        try:
            os.remove(temp_file)
        except Exception as e:
            logger.warning(f"Failed to delete temporary file {temp_file}: {e}")

    return starcode_exploded


def extract_summary(stdout: str) -> str:
    """ 
    Extracts the summary from the stdout of a bbduk2 command.
    
    Parameters:
        stdout (str): The stdout of the bbduk2 command
        
    Returns:
        str: The summary of the bbduk2 command
    """
    match = re.search(r"(Input:.*?Result:.*?bases\s*\(\d+\.\d+%\))", stdout, re.DOTALL)
    if match:
        return match.group(1)
    return None


def load_barcodes_chunked(barcodes_file: str, chunk_size: int):
    """
    Load barcodes in chunks from FASTQ files.

    Parameters:
        barcodes_file (str): Path to the barcodes FASTQ file.
        chunk_size (int): Number of records to read in each chunk.

    Yields:
        tuple: A chunk of fragment reads and barcode reads as lists of SeqRecords.
    """
    with gzip.open(barcodes_file, "rt") as bc_handle:
        bc_iter = SeqIO.parse(bc_handle, "fastq")
        while True:
            bc_chunk = list(islice(bc_iter, chunk_size))
            if not bc_chunk:
                break
            yield bc_chunk
            
            
def create_full_table(reads_BC: list)-> pd.DataFrame:
    """
    Create a full table of all found fragments with their corresponding barcodes and matching LUT data.
    
    Parameters:
        reads_BC (list): List of barcodes as SeqRecords
    
    Returns:
        pd.DataFrame: DataFrame containing the full table barcodes and ids
    """
    # Create full table with all found fragments and barcodes combinations
    full_table = pd.DataFrame({
        'ID': [record.id for record in reads_BC],
        'BC': [str(rec.seq) for rec in reads_BC]
    })

    return full_table


def analyze_tissue(file_path:str, data_dir:str, db:str, out_dir:str, library_fragments: pd.DataFrame,
                    lut_dna: pd.DataFrame, threads:int, bc_len:int, bbduk2_args: list, chunk_size:int) -> dict:
    """
    Analyze a single tissue sample based on its index in the load list.

    Parameters:
        file_path (str): The path to the FASTQ file
        data_dir (str): The directory containing the FASTQ files
        db (str): The path to the known barcodes FASTA file
        out_dir (str): The directory to save the found fragments
        library_fragments (pd.DataFrame): The library fragments
        lut_dna (pd.DataFrame): The LUT data
        threads (int): The number of threads to use
        bc_len (int): The length of the barcodes
        bbduk2_args (list): The arguments for bbduk2.sh
        
    Returns:
        dict: A log entry with the results
    """
    log_entry = {}
    
    # create the logname from the file_path
    log_entry['Name'] = os.path.basename(file_path).replace('.fastq.gz', '')
    
    path = os.path.join(data_dir, file_path)
    
    logger.info(f"Processing {file_path}")

    # Extraction of barcodes
    # ============================
    # Create a temporary file for the outputs
    out_name_BC = tempfile.NamedTemporaryFile(prefix="BC_", suffix=".fastq.gz", delete=False).name
    out_name_seqkit = tempfile.NamedTemporaryFile(prefix="BC_", suffix=".fastq.gz", delete=False).name
    
    # Run vsearch to extract the matching fragment reads fromt he tissue sample
    with tempfile.NamedTemporaryFile(delete=True, mode='w+', suffix='.txt') as vsearch_out, \
        tempfile.NamedTemporaryFile(delete=True, mode='w+', suffix='.txt') as keep_fragments:
        
        # Run vsearch and store the output in a temporary file allow for 1 missmatch
        vsearch_command = [
            f"zcat {path} | "
            f"vsearch --usearch_global - "
            f"--db {db} "
            "--id 0.95 "
            f"--blast6out {vsearch_out.name} "
            f"--threads {threads} "
            f"--minseqlength {bc_len}"
        ]
        _, stderr = run_command(vsearch_command, "vsearch", shell=True)

        match = re.search(r"(\d+ of \d+ \(\d+\.\d+%\))", stderr)
        if match:
            info = match.group(0)
            logger.info(f"Number of found Fragment reads that match to the reference: {info}")

        # Use awk to extract the matching fragment reads from the vsearch output
        awk_command = [f"awk '{{print $1}}' {vsearch_out.name} > {keep_fragments.name}"]
        _, stderr = run_command(awk_command, "awk", shell=True)
        
        # Extracting fragment reads and corresponding barcode reads that match to the reference

        # Use seqkit to extract the matching fragment reads
        seqkit_command = [
            f"seqkit grep -f {keep_fragments.name} {path} "
            f"-o filtered_fragments.fastq.gz -j {threads}"
        ]
        _, stderr = run_command(seqkit_command, "seqkit grep", shell=True)

        # Move the filtered fragment reads to the output directory
        shutil.move("filtered_fragments.fastq.gz", out_name_seqkit)

    # Run the bbduk2 command to extract barcodes from the barcode reads that match to a found fragment
    bbduk2_command = ["bbmap/bbduk2.sh"] + bbduk2_args + [
        f"threads={threads}",
        f"in={out_name_seqkit}",
        f"out={out_name_BC}",
    ]
    _ , stderr = run_command(bbduk2_command, f"bbduk2 barcode extraction for {file_path}")
    
    # Extract summary from the fragments extraction
    summary = extract_summary(stderr)
    if summary:
        logger.info(f"bbduk2 extraction summary:\n{summary}")

    # Save the Barcode extraction result in the log entry
    match = re.search(r"Result:\s*(\d+)", stderr)
    if match:
        log_entry['BC_reads'] = int(match.group(1))
    else:
        log_entry['BC_reads'] = 0

    # Chunked reading of barcodes
    # ============================
    # temp h5 output file_name
    output_file = os.path.join(out_dir, f"barcodes.h5")
    write_mode = 'w'
    
    for bc_chunk in load_barcodes_chunked(out_name_BC, chunk_size):
        # Create a full table of all found fragments with their corresponding barcodes and matching LUT data
        chunk_table = create_full_table(bc_chunk)
        
        # Save the full table to an h5 file
        chunk_table.to_hdf(output_file, key='data', mode=write_mode, format='table', append=(write_mode == 'a'))

        # After the first write, change mode to 'append'
        write_mode = 'a'
        # Explicitly free memory
        del bc_chunk, chunk_table
        gc.collect()
    
    # read in the full table
    try:
        barcode_table = pd.read_hdf(output_file, key='data')
        # remove the h5 file
        os.remove(output_file)
    except Exception as e:
        log_entry['unique_BCs'] = 0
        log_entry['scBCs'] = 0
        log_entry['SCdroppedBC'] = 0
        logger.info(f"Finished processing {file_path} found: {log_entry['BC_reads']} barcode reads; {log_entry['unique_BCs']} unique barcodes; {log_entry['scBCs']} reduced barcodes")
        return log_entry
    
    # Starcode based barcode reduction
    # ============================
    table_BC_sc = starcode_based_barcode_reduction(barcode_table)

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
    log_entry['unique_BCs'] = all_BCs
    log_entry['scBCs'] = sc_BCs
    log_entry['SCdroppedBC'] = SC_dropped_BC
    # Delete the old barcodes column
    barcode_table.drop(columns=['oldBC'], inplace=True)
    
    # Matching barcodes with information form the LUT and adding RNAcount
    # ============================
    # Reorganize the barcode_table to have BC and RNAcount
    BCcount = barcode_table['BC'].value_counts().reset_index()
    BCcount.columns = ['BC', 'RNAcount']
    # Extract only BC that are in BCcount
    foundFrags = library_fragments.merge(BCcount, on='BC', how='inner')
    if lut_dna is not None:
        # Merge with lut_dna on 'LUTnr'
        foundFrags = foundFrags.merge(lut_dna, on=['LUTnr','Peptide'], how='inner')
        # Rename the 'Reads' coulmn to 'Sequence'
        foundFrags.rename(columns={'Reads': 'Sequence'}, inplace=True)
    
    # Save the found fragments
    # ============================
    foundFrags.sort_values(by='RNAcount', ascending=False, inplace=True)
    output_filename = os.path.join(out_dir, f"found.{log_entry['Name']}.csv")
    
    foundFrags.to_csv(output_filename, index=False)
    
    logger.info(f"Finished processing {file_path} found: {log_entry['BC_reads']} barcode reads; {log_entry['unique_BCs']} unique barcodes; {log_entry['scBCs']} reduced barcodes")    
    
    return log_entry


def main():
    start_time = datetime.now()
    threads = multiprocessing.cpu_count()
    
    # load config
    config = get_config("S4")
    
    # Create a logger
    create_logger(config["log_dir"], "S4")

    # Try to load the necessary data
    try:
        library_fragments = pd.read_csv(config["input_table"], dtype={7: str})
        if config["in_name_LUT"] is not None:
            lut_dna = pd.read_csv(config["in_name_LUT"])
            # Drop the 'Sequence' column from lut_dna
            lut_dna.drop(columns='Sequence', inplace=True)
        else:
            lut_dna = None
    except Exception as e:
        logger.error(f"Error loading data: {e}")
        sys.exit(1)
    # Try to load the sample inputs
    try:
        load_list = pd.read_csv(config["sample_inputs"])
    except Exception as e:
        logger.error(f"Error loading sample inputs: {e}")
        sys.exit(1)
    # get the data directory and output directory from the config
    data_dir = config["sample_directory"]
    output_dir = config["output_dir"]
    db = config["db"]
    bc_len = config["bc_len"]
    # Create the output directory if it does not exist
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if not os.path.isdir(data_dir):
        logger.error(f"Data directory {data_dir} does not exist")
        sys.exit(1)
    log_table = []
    # get the settings for the barcode extraction
    bbduk2_args_BC = config["bbduk2_args"]
    chunk_size = config["chunk_size"]

    # Analyze each tissue sample
    for row in load_list.iterrows():
        # Extract the file name from the first column
        file_path = row[1]['Sample']
        log_entry = analyze_tissue(file_path, data_dir, db, output_dir, library_fragments, lut_dna, threads, bc_len, bbduk2_args_BC, chunk_size)
        if log_entry:
            log_table.append(log_entry)

    # Create a DataFrame from the log table
    log_df = pd.DataFrame(log_table)
    # Save the log table
    log_df.to_csv(config["log_file_path"], index=False)

    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
