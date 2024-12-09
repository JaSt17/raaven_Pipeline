#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script maps the found fragments and barcodes to the LUT data.
It then reduces the barcodes using the Starcode algorithm. And replaces the barcodes with the Starcode-reduced versions.
Finally it devides the barcodes into single-read and multi-read barcodes and splits the multi-read barcodes into clean and chimeric barcodes.
Saves all found barcodes to a CSV file.

Workflow:
    - Load the LUT data
    - Load the fragments and barcodes
    - Create a BLAST database from the LUT sequences
    - Save unique fragments to a FASTA file
    - Align unique fragments against the LUT database using BLASTn
    - Read the BLAST output into a DataFrame
    - Map every read to its corresponding LUTnr
    - Create a full table of all fragments that matched the LUT with their BLASTn results
    - Perform barcode reduction using Starcode clustering
    - Replace barcodes with Starcode-reduced versions
    - Split reads into single-read and multi-read barcodes
    - Split multi-read barcodes into clean and chimeric ones. (chimeric barcodes means one barcode maps to multiple LUTnr)
    - Create consensus alignment of chimeric barcodes (get the barcode with the highest maximal read count for each LUTnr)
    - Combine all tables into final output

Inputs for the script are:
    - in_name_LUT: Path to the LUT file
    - fragment_file: Path to the fragments FASTQ file
    - barcode_file: Path to the barcodes FASTQ file
    - out_name: Path to the output file

Output of the script is:
    - A CSV file containing the found barcodes with the following columns:
        BC,LUTnr,bitScore,tCount,mCount,Mode
"""

import gzip
import os
import tempfile
import subprocess
import multiprocessing
import pandas as pd
from Bio import SeqIO
from itertools import islice
from Bio.Seq import translate
from datetime import datetime
import logging
import sys
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


def run_command(command: list, description: str, shell=False) -> tuple:
    """
    Runs a subprocess command and returns stdout, stderr, and error status.
    
    Parameters: 
        command (list or str): The command to run
        description (str): A description of the command
        shell (bool): Whether to execute the command through the shell
    
    Returns:   
        tuple: stdout and stderr
    """
    # change the command to a string if it is a list if shell is True
    if shell:
        command = command[0]
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=shell) as process:
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            logger.error(f"Error running {description} with code {process.returncode}")
            logger.error(stderr)
            sys.exit(1)
        return stdout, stderr
    
    
def save_unique_fragments_barcodes(fragments_file: str, barcodes_file) -> tuple:
    """
    Save all found unique fragments and barcodes from the library to a FASTA file.
    
    Parameters:
        reads_frag (list): List of fragments as SeqRecords
        
    Returns:
        tuple: Path to the FASTA file containing unique fragments and the unique fragments
    """
    logger.info("Saving unique fragments and barcodes to FASTA file")

    # Get the number of threads
    num_threads = multiprocessing.cpu_count()
    
    out_name_1 = "/".join(fragments_file.split("/")[:-1]) + "/unique_fragments.fasta"
    # Build shell command for extracting unique sequences
    command = [
        f"zcat {fragments_file} | "
        "awk 'NR % 4 == 2' | "
        f"sort --parallel={num_threads} --buffer-size=1G | "
        f"uniq | grep -v 'N' | "
        "awk '{print \">\" NR \"\\n\" $0}' > "
        f"{out_name_1}"
    ]
    # Execute shell command
    run_command(command, "Extract unique sequences", shell=True)
    
    number_of_unique_fragments, _ = run_command([f" echo $(( $(wc -l < {out_name_1}) / 2 ))"], "Extract unique sequences", shell=True)
    logger.info(f"Number of unique fragments: {number_of_unique_fragments.strip()}")
    
    out_name_2 = "/".join(barcodes_file.split("/")[:-1]) + "/unique_barcodes.fasta"
    # Build shell command for extracting unique sequences
    command = [
        f"zcat {barcodes_file} | "
        "awk 'NR % 4 == 2' | "
        f"sort --parallel={num_threads} --buffer-size=1G | "
        f"uniq | grep -v 'N' | "
        "awk '{print \">\" NR \"\\n\" $0}' > "
        f"{out_name_2}"
    ]
    # Execute shell command
    run_command(command, "Extract unique sequences", shell=True)
    
    number_of_unique_barcodes, _ = run_command([f" echo $(( $(wc -l < {out_name_2}) / 2 ))"], "Extract unique sequences", shell=True)
    logger.info(f"Number of unique barcodes: {number_of_unique_barcodes.strip()}")

    return out_name_1, out_name_2


def load_frag_bc_reads_chunked(fragments_file: str, barcodes_file: str, chunk_size: int):
    """
    Load fragments and barcodes in chunks from FASTQ files.

    Parameters:
        fragments_file (str): Path to the fragments FASTQ file.
        barcodes_file (str): Path to the barcodes FASTQ file.
        chunk_size (int): Number of records to read in each chunk.

    Yields:
        tuple: A chunk of fragment reads and barcode reads as lists of SeqRecords.
    """
    with gzip.open(fragments_file, "rt") as frag_handle, gzip.open(barcodes_file, "rt") as bc_handle:
        frag_iter = SeqIO.parse(frag_handle, "fastq")
        bc_iter = SeqIO.parse(bc_handle, "fastq")
        while True:
            frag_chunk = list(islice(frag_iter, chunk_size))
            bc_chunk = list(islice(bc_iter, chunk_size))
            if not frag_chunk or not bc_chunk:
                break
            yield frag_chunk, bc_chunk


def create_full_table(reads_frag: list, reads_BC: list)-> pd.DataFrame:
    """
    Create a full table with fragments and barcodes.
    
    Parameters:
        reads_frag (list): List of fragments as SeqRecords
        reads_BC (list): List of barcodes as SeqRecords
    
    Returns:
        pd.DataFrame: DataFrame containing the full table with fragments and barcodes
    """
    # Create full table with fragments and barcodes pair
    full_table = pd.DataFrame({
        'Reads': [str(rec.seq) for rec in reads_frag],
        'BC': [str(rec.seq) for rec in reads_BC]
    })
    # create unique Lutnrs for each fragment and translate the fragment to a peptide
    # create a DataFrame thatonly contains unique fragments
    unique_fragments = full_table['Reads'].unique()
    lutnrs_df = pd.DataFrame({'Reads': unique_fragments})
    # add the Peptide column to the DataFrame
    lutnrs_df['Peptide'] = lutnrs_df['Reads'].apply(lambda x: translate(x[2:-2]))
    # change all * to M in the Peptide column
    lutnrs_df['Peptide'] = lutnrs_df['Peptide'].str.replace('*', 'M')
    # add a column with the LUTnr for each fragment
    lutnrs_df['LUTnr'] = lutnrs_df.index
    # add seq_ in front of the LUTnr
    lutnrs_df['LUTnr'] = 'seq_' + lutnrs_df['LUTnr'].astype(str)
    
    # merge the lutnrs_df with the full table
    full_table = full_table.merge(lutnrs_df, on='Reads', how='left')

    return full_table


def starcode_based_reduction_and_replace(full_table: pd.DataFrame, input_file_name: str, columns_name: str)-> pd.DataFrame:
    """
    Perform reduction using Starcode clustering and replaces with the reduced versions.

    Parameters:
        full_table (pd.DataFrame): DataFrame containing the full table with fragments and barcodes
        columns_name (str): Name of the column containing the sequences
        
    Returns:
        pd.DataFrame: DataFrame containing the Starcode-reduced barcodes
    """
    logger.info("Running Starcode clustering")
    
    # get the number of threads
    num_threads = multiprocessing.cpu_count()
    starcode_output = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".txt")
    starcode_cmd = f"gunzip -c {input_file_name} | starcode -t {num_threads-1} --print-clusters -d 1 -r5 -q -o {starcode_output.name}"
    subprocess.run(starcode_cmd, shell=True, check=True)

    starcode_column = f"sc{columns_name}"
    # Read Starcode output
    starcode_df = pd.read_csv(starcode_output.name, sep='\t', header=None, names=[starcode_column, 'count', 'seq_list'])
    starcode_df['seq_list'] = starcode_df['seq_list'].str.split(',')
    number_of_clusters = len(starcode_df)

    # Explode the seq_list to have one row per seq and rename the columns
    starcode_exploded = starcode_df.explode('seq_list')
    starcode_exploded.rename(columns={'seq_list': columns_name}, inplace=True)
    # only keep the unique barcodes
    starcode_exploded.drop_duplicates(subset=[columns_name], inplace=True)
    logger.info(f"Number of unique {columns_name} after Starcode reduction: {number_of_clusters}")
    
    # combine the full table with the starcode exploded table based on the input column
    logger.info(f"Replacing {columns_name} with Starcode-reduced versions")
    full_table = full_table.merge(starcode_exploded[[columns_name, starcode_column]], on=columns_name, how='left')
    # rename the coulmns and drop the old BC column
    full_table.rename(columns={columns_name: 'old', starcode_column: columns_name}, inplace=True)
    full_table.drop(columns=['old'], inplace=True)
    
    return full_table


def split_reads_into_single_and_multi_read_barcodes(full_table: pd.DataFrame)-> tuple:
    """
    Split reads into single-read and multi-read barcodes.

    Parameters:
        full_table (pd.DataFrame): DataFrame containing the full table with fragments and barcodes
        
    Returns:
        tuple: DataFrames containing the single-read and multi-read barcodes
    """
    logger.info("Splitting reads into single-read and multi-read barcodes")
    barcode_counts = full_table['BC'].value_counts()
    # get single read barcodes table
    single_read_barcodes = barcode_counts[barcode_counts == 1].index
    # get multi read barcodes table
    multi_read_barcodes = barcode_counts[barcode_counts > 1].index

    logger.info(f"Number of single-read barcodes: {len(single_read_barcodes)}")
    logger.info(f"Number of multi-read barcodes: {len(multi_read_barcodes)}")

    # create tables with single and multi read barcodes
    temp_table_single = full_table[full_table['BC'].isin(single_read_barcodes)].copy()
    temp_table_multi = full_table[full_table['BC'].isin(multi_read_barcodes)].copy()

    # set the max read count, total read count, and mode for single read barcodes
    temp_table_single['mCount'] = 1
    temp_table_single['tCount'] = 1
    temp_table_single['Mode'] = 'Amb'

    return temp_table_single, temp_table_multi


def split_multi_read_barcodes_into_clean_and_chimeric(temp_table_multi: pd.DataFrame)-> tuple:
    """
    Split multi-read barcodes into clean and chimeric based on LUTnr.

    Parameters:
        temp_table_multi (pd.DataFrame): DataFrame containing the multi-read barcodes
        
    Returns:
        tuple: DataFrames containing the clean and chimeric multi-read barcodes
    """
    logger.info("Splitting multi-read barcodes into clean and chimeric")

    # Group by BC and LUTnr to compute statistics for every single pair
    temp_table_multi_grouped = temp_table_multi.groupby(['BC', 'LUTnr']).agg({
        'Reads': 'count'
    }).reset_index().rename(columns={'Reads': 'tCount'}) # tCount is the total count of reads for a given BC and LUTnr pair

    # Identify clean and chimeric barcodes
    bc_counts = temp_table_multi_grouped['BC'].value_counts()
    # if the barcode only appears once, then it is a clean barcode since it clearly maps to one LUTnr
    clean_barcodes = bc_counts[bc_counts == 1].index
    # if the barcode appears more than once, then it is a chimeric barcode since it maps to multiple LUTnr
    chimeric_barcodes = bc_counts[bc_counts > 1].index

    # create tables with clean and chimeric barcodes
    temp_table_multi_clean = temp_table_multi_grouped[temp_table_multi_grouped['BC'].isin(clean_barcodes)].copy()
    temp_table_multi_chimeric = temp_table_multi_grouped[temp_table_multi_grouped['BC'].isin(chimeric_barcodes)].copy()

    # set the max read count, total read count, since we only have one barcode for this LUTnr and mode for clean barcode
    temp_table_multi_clean['mCount'] = temp_table_multi_clean['tCount']
    temp_table_multi_clean['Mode'] = 'Def'

    return temp_table_multi_clean, temp_table_multi_chimeric


def get_valid_chimeric_barcodes(temp_table_multi_chimeric: pd.DataFrame, threshold: float)-> pd.DataFrame:
    """
    Calculate consensus alignment of chimeric barcodes. This means we are getting the barcode with the highest maximal read count for each LUTnr.

    Parameters:
        temp_table_multi_chimeric (pd.DataFrame): DataFrame containing the chimeric multi-read barcodes
    
    Returns:
        pd.DataFrame: DataFrame containing the consensus alignment for chimeric barcodes
    """
    logger.info(f"Extracting chimeric barcodes with maximal read ratio above {threshold}")
    # For every barcode LUTnr pair, set mCount to tCount
    temp_table_multi_chimeric['mCount'] = temp_table_multi_chimeric['tCount']
    # Sum the total count of reads for each barcode
    temp_table_multi_chimeric['tCount'] = temp_table_multi_chimeric.groupby('BC')['tCount'].transform('sum')
    temp_table_multi_chimeric['Mode'] = 'Amb'
    
    idx = temp_table_multi_chimeric.groupby('BC')['mCount'].idxmax()
    temp_table = temp_table_multi_chimeric.loc[idx].copy()
    
    # get the index of all rows where the maximal read count divided by the total read count is above the threshold
    idx = temp_table['mCount'] / temp_table['tCount'] >= threshold
    
    # set the mode to Def for all rows where the maximal read count divided by the total read count is above the threshold
    temp_table.loc[idx, 'Mode'] = 'Def'

    return temp_table


def combine_tables(temp_table_multi_clean: pd.DataFrame, temp_table_multi_consensus: pd.DataFrame, temp_table_single: pd.DataFrame, threshold:float)-> pd.DataFrame:
    """
    Combine the multi-read and single-read tables into the final output table.
    
    Parameters:
        temp_table_multi_clean (pd.DataFrame): DataFrame containing the clean multi-read barcodes
        temp_table_multi_consensus (pd.DataFrame): DataFrame containing the consensus alignment for chimeric barcodes
        temp_table_single (pd.DataFrame): DataFrame containing the single-read barcodes
        
    Returns:
        pd.DataFrame: DataFrame containing the final output table
    """
    logger.info("Combining tables to create final output")
    # Combine clean and consensus tables
    temp_table_multi_final = pd.concat([temp_table_multi_clean, temp_table_multi_consensus], ignore_index=True)
    logger.info(f"Number of barcodes-fragment pairs sequenced more than once: {len(temp_table_multi_final)}")
    logger.info(f"  Number of barcodes mapping to only one fragment: {len(temp_table_multi_clean)}")
    logger.info(f"  Number of barcodes mapping to more than one fragment: {len(temp_table_multi_consensus)}")
    logger.info(f"      From these {len(temp_table_multi_consensus[temp_table_multi_consensus['Mode'] == 'Def'])} are clean barcodes (ratio above {threshold})")
    logger.info(f"      From these {len(temp_table_multi_consensus[temp_table_multi_consensus['Mode'] == 'Amb'])} are chimeric barcodes (ratio below {threshold})")
    logger.info(f"Number of barcodes-fragment pairs sequenced only once: {len(temp_table_single)}")
    output_table = pd.concat([temp_table_multi_final, temp_table_single], ignore_index=True)

    return output_table


def main():
    start_time = datetime.now()

    # load config
    config = get_config("S3")
    
    # Create a logger
    create_logger(config["log_dir"], "S3")
    
    # Temporary directory for intermediate results
    temp_dir = tempfile.mkdtemp()
    chunk_files = []

    # Chunk size for reading fragments and barcodes
    chunk_size = config["chunk_size"]

    # Process fragments and barcodes in chunks
    for i, (frag_chunk, bc_chunk) in enumerate(load_frag_bc_reads_chunked(config["fragment_file"], config["barcode_file"], chunk_size)):
        logger.info(f"Processing chunk {i + 1}")
        
        # Create a full table with fragments and barcodes from the chunk
        chunk_table = create_full_table(frag_chunk, bc_chunk)

        # Save chunk results to temporary file
        chunk_file = os.path.join(temp_dir, f"chunk_{i}.pkl")
        chunk_table.to_pickle(chunk_file)
        chunk_files.append(chunk_file)

    # Combine all chunk tables
    logger.info("Combining all chunks into the full table")
    full_table = pd.concat([pd.read_pickle(chunk_file) for chunk_file in chunk_files], ignore_index=True)
    del chunk_files
    
    save_unique_fragments_barcodes(config["fragment_file"], config["barcode_file"])
    
    # Perform fragment reduction using Starcode clustering
    full_table = starcode_based_reduction_and_replace(full_table, config['fragment_file'], 'Reads')
    
    # Perform barcode reduction using Starcode clustering
    full_table = starcode_based_reduction_and_replace(full_table, config['barcode_file'], 'BC')

    # Split reads into single-read and multi-read barcodes
    temp_table_single, temp_table_multi = split_reads_into_single_and_multi_read_barcodes(full_table)
    del full_table

    # Split multi-read barcodes into clean and chimeric
    temp_table_multi_clean, temp_table_multi_chimeric = split_multi_read_barcodes_into_clean_and_chimeric(temp_table_multi)

    # Create consensus alignment for chimeric barcodes and get valid chimeric barcodes based on threshold
    temp_table_multi_chimeric = get_valid_chimeric_barcodes(temp_table_multi_chimeric, config["threshold"])

    # Combine all tables into final output
    output_table = combine_tables(temp_table_multi_clean, temp_table_multi_chimeric, temp_table_single, config["threshold"])
    del temp_table_multi_clean, temp_table_multi_chimeric, temp_table_single

    # change reads to Sequence
    output_table.rename(columns={'Reads': 'Sequence'}, inplace=True)

    # Save the output table
    output_table.to_csv(config['out_name'], index=False)
    logger.info(f"Output saved to {config['out_name']}")

    # Print total analysis time
    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
