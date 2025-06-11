#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script extracts barcodes from given samples, and matches them with corresponding fragments suing the LUT.
It processes a csv file with Sample and Group and saves the results in a log table and saves the found fragments in a csv file.

Workflow:
    - Load the library fragments and LUT data
    - For each RNA sample:
        - Extract barcodes using bbduk2.sh
        - Match barcodes with reference barcodes using vsearch
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
import logging
from datetime import datetime
import pandas as pd
import multiprocessing
from pathlib import Path
import gzip
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
        bc_iter = SeqIO.parse(bc_handle, "fasta")
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


def write_fasta(df, path, id_col='BC', seq_col='BC'):
    """    Write a DataFrame to a FASTA file.
    Parameters:
        df (pd.DataFrame): DataFrame containing the sequences
        path (str): Path to save the FASTA file
        id_col (str): Column name for sequence IDs
        seq_col (str): Column name for sequences
    """
    with open(path, 'w') as f:
        for i, row in df.iterrows():
            f.write(f'>{row[id_col]}\n{row[seq_col]}\n')


def parse_blast6(blast_path):
    """ Parse a BLAST6 output file into a DataFrame.
    Parameters:
        blast_path (str): Path to the BLAST6 output file
    """
    cols = [
        'query_BC', 'matched_BC', 'pident', 'length',
        'mismatch', 'gapopen', 'qstart', 'qend',
        'sstart', 'send', 'evalue', 'bitscore'
    ]
    df = pd.read_csv(blast_path, sep='\t', names=cols)
    return df

def starcode_merge(unique_barcodes, barcode_db, threshold=0.95):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)

        # Paths
        query_fasta = tmpdir / 'queries.fasta'
        db_fasta = tmpdir / 'database.fasta'
        blast_out = tmpdir / 'matches.b6'

        # Write FASTA files
        write_fasta(unique_barcodes, query_fasta)
        write_fasta(barcode_db, db_fasta)
        
        # get the number of threads
        threads = os.cpu_count()

        # Run vsearch
        cmd = [
            'vsearch',
            '--usearch_global', str(query_fasta),
            '--db', str(db_fasta),
            '--id', str(threshold),
            '--blast6out', str(blast_out),
            '--threads', f'{threads}',
        ]
        subprocess.run(cmd, check=True)

        # Parse results
        result_df = parse_blast6(blast_out)
        # rename 'matched_BC' to 'BC' for consistency
        result_df.rename(columns={'matched_BC': 'BC'}, inplace=True)
        
        found_barcodes = unique_barcodes.merge(result_df, on='BC', how='inner')
        
        return found_barcodes


def analyze_tissue(file_path:str, data_dir:str, db:str, starcode:bool, out_dir:str, library_fragments: pd.DataFrame,
                    lut_dna: pd.DataFrame, threads:int, bbduk2_args: list,) -> dict:
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
    
    # Check if the output file already exists
    out_file = os.path.join(out_dir, f"unique_barcodes_{log_entry['Name']}.csv")
    if os.path.isfile(out_file):
        logger.info(f"Output file {out_file} already exists.")
        unique_barcodes = pd.read_csv(out_file)
        if not os.path.isfile(db):
            logger.error(f"Database file {db} does not exist.")
            sys.exit(1)
        # read in all barcode form the db.fasta
        barcode_dict = SeqIO.to_dict(SeqIO.parse(db, "fasta"))
            # create a DataFrame from the barcode_db
        barcode_db = pd.DataFrame({
            'BC': [str(record.seq) for record in barcode_dict.values()],
        })

    else:
        # Extraction of barcodes first
        # ============================
        # Create a temporary file for the barcode outputs
        out_name_BC = tempfile.NamedTemporaryFile(prefix="BC_", suffix=".fastq.gz", delete=False).name

        # Run the bbduk2 command to extract barcodes first
        bbduk2_command = ["bbmap/bbduk2.sh"] + bbduk2_args + [
            f"threads={threads}",
            f"in={path}",
            f"out={out_name_BC}",
        ]
        _, stderr = run_command(bbduk2_command, f"bbduk2 barcode extraction for {path}")

        # Extract summary from the barcode extraction
        summary = extract_summary(stderr)
        if summary:
            logger.info(f"bbduk2 extraction summary:\n{summary}")
        total_reads = int(re.search(r"Input:\s+(\d+)", summary).group(1))
        per_mil_scale = total_reads / 1_000_000
            
        # Extract the the number of unique barcodes found
        stdout, _ = run_command([f"zcat {out_name_BC} | awk 'NR % 4 == 2' | sort | uniq -c | sort -nr "], "Count unique barcodes", shell=True)
        # covert the stdout to an table with count and barcode
        unique_barcodes = pd.DataFrame([line.split() for line in stdout.strip().split('\n')], columns=['Count', 'BC'])
        # scale the counts to per million reads
        unique_barcodes['Count'] = unique_barcodes['Count'].astype(int)
        unique_barcodes['Count_per_mil_reads'] = unique_barcodes['Count'].astype(int) / per_mil_scale
        # save all found unique barcodes in a csv file with their counts
        unique_barcodes.to_csv(os.path.join(out_dir, f"unique_barcodes_{log_entry['Name']}.csv"), index=False)
        logger.info(f"Number of unique barcodes found: {unique_barcodes.shape[0]}")
        
        # match the barcodes with the reference barcodes from the plasmid database
        if not os.path.isfile(db):
            logger.error(f"Database file {db} does not exist.")
            sys.exit(1)
        # read in all barcode form the db.fasta
        barcode_dict = SeqIO.to_dict(SeqIO.parse(db, "fasta"))
            # create a DataFrame from the barcode_db
        barcode_db = pd.DataFrame({
            'BC': [str(record.seq) for record in barcode_dict.values()],
        })

    try:
        # extract information about barcodes found in the sampes
        log_entry['BC_reads'] = int(unique_barcodes['Count'].sum())
        log_entry['unique_BC'] = int(unique_barcodes.shape[0])
        
        # Merge the found barcodes with the barcode database
        if starcode:
            # Use starcode to merge similar barcodes
            BCcount = starcode_merge(unique_barcodes, barcode_db)
        else:
            BCcount = unique_barcodes.merge(barcode_db, on='BC', how='inner')
            
        # extract information about matched barcodes
        log_entry['matched_BC_reads'] = int(BCcount['Count'].sum())
        log_entry['matched_unique_BC'] = int(BCcount.shape[0])
        
        #rename the columns to 'BC' and 'Count'
        BCcount.rename(columns={'BC': 'BC', 'Count': 'RNAcount', 'Count_per_mil_reads': 'RNAcount_per_mil_reads'}, inplace=True)
            
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
        
        # Log the information for each sample
        logger.info(
            f"Finished processing {file_path} found: "
            f"{log_entry['BC_reads']} barcode reads; "
            f"{log_entry['unique_BC']} unique barcodes; "
            f"{log_entry['matched_BC_reads']} matched barcode reads; "
            f"{log_entry['matched_unique_BC']} matched unique barcodes; ")
            
        return log_entry
    
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        logger.error(f"Error processing {file_path}: {e}")
        exit(1)
    
def create_summary_csv(directory):
    """
    Create a summary CSV from all unique barcode files in the specified directory.
    
    Parameters:
    - directory: str, path to the directory containing the unique barcode CSV files.
    
    Returns:
    - None, saves the summary CSV to the specified directory.
    """
    
    # Ensure the directory exists
    if not os.path.exists(directory):
            raise FileNotFoundError(f"The directory {directory} does not exist.")
        
    # List to hold data from all files
    all_barcodes = []

    # Loop through all relevant CSV files
    for filename in os.listdir(directory):
        if filename.startswith('unique_barcodes') and filename.endswith('.csv'):
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path)

            # Extract sample name from filename
            sample_name = filename.split('unique_barcodes_')[1].split('.csv')[0]

            # Filter for Count > 1000
            df = df[df['Count'] > 100].copy()

            # Add sample name column
            df['Sample'] = sample_name

            all_barcodes.append(df)

    # Concatenate all sample DataFrames
    all_barcodes_df = pd.concat(all_barcodes, ignore_index=True)

    # --- Pivot the data to create one column per sample ---
    pivot_df = all_barcodes_df.pivot_table(
        index='BC',
        columns='Sample',
        values='Count',
        aggfunc='sum',
        fill_value=0
    )

    # --- Add a total 'Count' column (sum across samples) ---
    pivot_df['Count'] = pivot_df.sum(axis=1)

    # --- Add 'Samples_Found_In' column: number of samples with a count > 0 ---
    sample_columns = [col for col in pivot_df.columns if col != 'Count']
    pivot_df['Samples_Found_In'] = pivot_df[sample_columns].gt(0).sum(axis=1)

    # Reorder columns: Count, Samples_Found_In, then sample columns
    cols = ['Count', 'Samples_Found_In'] + sample_columns
    pivot_df = pivot_df[cols]

    # Sort by total Count descending
    pivot_df = pivot_df.sort_values(by=['Samples_Found_In','Count'], ascending=[False, False])

    # Save to CSV
    pivot_df.to_csv(f'{directory}.csv')


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

    # Analyze each tissue sample
    for row in load_list.iterrows():
        # Extract the file name from the first column
        file_path = row[1]['Sample']
        log_entry = analyze_tissue(file_path, data_dir, db, config["starcode"],
                                    output_dir, library_fragments, lut_dna, threads,
                                    bbduk2_args_BC)
        if log_entry:
            log_table.append(log_entry)

    # Create a DataFrame from the log table
    log_df = pd.DataFrame(log_table)
    # Save the log table
    log_df.to_csv(config["log_file_path"], index=False)
    
    # Create a summary CSV from all unique barcode files in the output directory
    create_summary_csv(config["output_dir"])

    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
