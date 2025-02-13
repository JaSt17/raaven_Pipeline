#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script reduces the barcodes using the Starcode. And replaces the barcodes with the Starcode-reduced versions.
Finally it devides the barcodes into single-read and multi-read barcodes and splits the multi-read barcodes into clean and chimeric barcodes.
Then we set all single-read barcodes to Mode 'Single' and all clean multi-read barcodes to Mode 'Def'.
The chimeric barcodes are set to Mode 'Chimeric' if the maximal read count divided by the total read count is below a certain threshold.
Saves all found barcodes to 3 different CSV files:
    - Definitiv barcodes
    - Chimeric barcodes
    - Single barcodes

Workflow:
    - Save unique fragments and barcodes to a FASTA file
    - Create a full table of all fragments that matched the LUT with their BLASTn results
        - This is done for every chunk of fragments and barcodes since the files are too large to process at once
    - Perform barcode reduction using Starcode clustering
    - Replace barcodes with Starcode-reduced versions
    - Split reads into single-read and multi-read barcodes
    - Split multi-read barcodes into clean and chimeric ones. (chimeric barcodes means one barcode maps to multiple LUTnr)
    - Set the chimeric barcodes to Mode 'Amb' if the maximal read count divided by the total read count is below a certain threshold
    - Combine all tables into final output

Inputs for the script are:
    - in_name_LUT: Path to the LUT file
    - fragment_file: Path to the fragments FASTQ file
    - barcode_file: Path to the barcodes FASTQ file
    - threshold: Threshold for the ratio of the most frequent barcode to all found barcodes for chimeric barcode detection
    - chunk_size: Number of records to read in each chunk
    - out_name: Path to the output file

Output of the script is:
    - A CSV file containing the found barcodes with the following columns:
        BC,LUTnr,tCount,Reads,Peptide,mCount,Mode
"""

import gzip
import gc
import os
import tempfile
import subprocess
import multiprocessing
from itertools import islice
import pandas as pd
from Bio import SeqIO
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


def save_unique_fragments_barcodes(fragments_file: str, barcodes_file: str) -> None:
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
        f"sort --parallel={num_threads}| "
        f"uniq | grep -v 'N' | "
        "awk '{print \">\" NR \"\\n\" $0}' > "
        f"{out_name_1}"
    ]
    # Execute shell command
    run_command(command, "Extract unique sequences", shell=True)
    
    number_of_unique_fragments, _ = run_command(
        [f"echo $(( $(cat {out_name_1} | wc -l) / 2 ))"], 
        "Extract unique sequences", 
        shell=True
    )
    logger.info(f"Number of unique fragment reads: {number_of_unique_fragments.strip()}")
    
    out_name_2 = "/".join(barcodes_file.split("/")[:-1]) + "/unique_barcodes.fasta"
    # Build shell command for extracting unique sequences
    command = [
        f"zcat {barcodes_file} | "
        "awk 'NR % 4 == 2' | "
        f"sort --parallel={num_threads}| "
        f"uniq | grep -v 'N' | "
        "awk '{print \">\" NR \"\\n\" $0}' > "
        f"{out_name_2}"
    ]
    # Execute shell command
    run_command(command, "Extract unique sequences", shell=True)
    
    number_of_unique_barcodes, _ = run_command(
        [f"echo $(( $(cat {out_name_2} | wc -l) / 2 ))"], 
        "Extract unique sequences", 
        shell=True
    )
    logger.info(f"Number of unique barcode reads: {number_of_unique_barcodes.strip()}")


def create_full_table(lut_df: pd.DataFrame, reads_frag: list, reads_BC: list)-> pd.DataFrame:
    """
    Create a full table of all found fragments with their corresponding barcodes and matching LUT data.
    
    Parameters:
        lut_df (pd.DataFrame): DataFrame containing the LUT sequences
        reads_frag (list): List of fragments as SeqRecords
        reads_BC (list): List of barcodes as SeqRecords
    
    Returns:
        pd.DataFrame: DataFrame containing the full table with fragments and barcodes
    """
    # Create full table with all found fragments and barcodes combinations
    full_table = pd.DataFrame({
        'Sequence': [str(rec.seq) for rec in reads_frag],
        'BC': [str(rec.seq) for rec in reads_BC]
    })
    # extract the Sequnce without linkers from the LUT
    lut_df['Sequence'] = lut_df['Sequence'].str.extract(r'([ATGC]+)')
    
    # Merge the full table with the LUT data based on the Sequence
    full_table = full_table.merge(lut_df[['Sequence', 'LUTnr', 'Peptide']], on='Sequence', how='inner')

    # rename the Fragment column to Reads
    full_table.rename(columns={'Sequence': 'Reads'}, inplace=True)

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
    
    total_number_of_barcodes = len(barcode_counts)

    logger.info(f"Number of single-read barcodes: {len(single_read_barcodes)} ({len(single_read_barcodes)/total_number_of_barcodes*100:.2f}%)")
    logger.info(f"Number of multi-read barcodes: {len(multi_read_barcodes)} ({len(multi_read_barcodes)/total_number_of_barcodes*100:.2f}%)")

    # create tables with single and multi read barcodes
    temp_table_single = full_table[full_table['BC'].isin(single_read_barcodes)].copy()
    temp_table_multi = full_table[full_table['BC'].isin(multi_read_barcodes)].copy()

    # set the max read count, total read count, and mode for single read barcodes
    temp_table_single['mCount'] = 1
    temp_table_single['tCount'] = 1
    temp_table_single['Mode'] = 'Single'

    return temp_table_single, temp_table_multi


def split_multi_read_barcodes_into_clean_and_chimeric(temp_table_multi: pd.DataFrame) -> tuple:
    """
    Split multi-read barcodes into clean and chimeric based on LUTnr, including Peptide information.

    Parameters:
        temp_table_multi (pd.DataFrame): DataFrame containing the multi-read barcodes
        
    Returns:
        tuple: DataFrames containing the clean and chimeric multi-read barcodes
    """
    logger.info("Splitting multi-read barcodes into clean and chimeric")

    # Group by BC and LUTnr to compute statistics and collect reads and peptides for every single pair
    temp_table_multi_grouped = temp_table_multi.groupby(['BC', 'LUTnr']).agg({
        'Reads': lambda x: list(x),  # Collect all reads in a list
        'Reads': 'count',           # Count the number of reads
        'Peptide': lambda x: list(x)  # Collect all peptides in a list
    }).reset_index().rename(columns={'Reads': 'tCount'})  # tCount is the total count of reads for a given BC and LUTnr pair

    # Identify clean and chimeric barcodes
    bc_counts = temp_table_multi_grouped['BC'].value_counts()
    # If the barcode only appears once, it is clean since it maps to one LUTnr
    clean_barcodes = bc_counts[bc_counts == 1].index
    # If the barcode appears more than once, it is chimeric since it maps to multiple LUTnr
    chimeric_barcodes = bc_counts[bc_counts > 1].index

    # Separate tables with clean and chimeric barcodes
    temp_table_multi_clean = temp_table_multi_grouped[temp_table_multi_grouped['BC'].isin(clean_barcodes)].copy()
    temp_table_multi_chimeric = temp_table_multi_grouped[temp_table_multi_grouped['BC'].isin(chimeric_barcodes)].copy()

    # Add the `Reads` and `Peptide` columns to clean and chimeric tables
    clean_reads = temp_table_multi[temp_table_multi['BC'].isin(clean_barcodes)].groupby(['BC', 'LUTnr']).agg({
        'Reads': lambda x: list(x),
        'Peptide': lambda x: list(x)
    }).reset_index()

    chimeric_reads = temp_table_multi[temp_table_multi['BC'].isin(chimeric_barcodes)].groupby(['BC', 'LUTnr']).agg({
        'Reads': lambda x: list(x),
        'Peptide': lambda x: list(x)
    }).reset_index()

    # Ensure the `Reads` and `Peptide` columns are unique lists for clean and chimeric barcodes
    clean_reads['Reads'] = clean_reads['Reads'].apply(lambda x: x[0])
    clean_reads['Peptide'] = clean_reads['Peptide'].apply(lambda x: x[0])

    chimeric_reads['Reads'] = chimeric_reads['Reads'].apply(lambda x: x[0])
    chimeric_reads['Peptide'] = chimeric_reads['Peptide'].apply(lambda x: x[0])
    
    # remove the 'Peptide' column from the tables before merging
    temp_table_multi_clean.drop(columns=['Peptide'], inplace=True)
    temp_table_multi_chimeric.drop(columns=['Peptide'], inplace=True)

    temp_table_multi_clean = pd.merge(temp_table_multi_clean, clean_reads, on=['BC', 'LUTnr'])
    temp_table_multi_chimeric = pd.merge(temp_table_multi_chimeric, chimeric_reads, on=['BC', 'LUTnr'])

    # Set additional metadata for clean barcodes
    temp_table_multi_clean['mCount'] = temp_table_multi_clean['tCount']  # Maximum count is total count for clean
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
    
    # sort the table by BC and mCount
    temp_table_multi_chimeric = temp_table_multi_chimeric.sort_values(by=['BC', 'mCount'], ascending=[True, False])
    
    temp_table_tcount = temp_table_multi_chimeric.groupby('BC').agg({
        'tCount': 'sum',
    }).reset_index()
    # remove the tCount column from the table
    temp_table_multi_chimeric.drop(columns=['tCount'], inplace=True)
    
    # merge the tCount column to the temp_table_multi_chimeric
    temp_table_multi_chimeric = pd.merge(temp_table_multi_chimeric, temp_table_tcount, on='BC')
    
    # Set the mode based on the threshold for the maximal read count divided by the total read count
    temp_table_multi_chimeric['Mode'] = temp_table_multi_chimeric.apply(lambda x: 'Def' if x['mCount'] / x['tCount'] >= threshold else 'Chimeric', axis=1)

    return temp_table_multi_chimeric


def combine_tables(temp_table_multi_clean: pd.DataFrame, temp_table_multi_chimeric: pd.DataFrame, temp_table_single: pd.DataFrame, threshold:float)-> pd.DataFrame:
    """
    Combine the multi-read and single-read tables into the final output table.
    
    Parameters:
        temp_table_multi_clean (pd.DataFrame): DataFrame containing the clean multi-read barcodes
        temp_table_multi_chimeric (pd.DataFrame): DataFrame containing the consensus alignment for chimeric barcodes
        temp_table_single (pd.DataFrame): DataFrame containing the single-read barcodes
        
    Returns:
        pd.DataFrame: DataFrame containing the final output table
    """
    # Combine clean and consensus tables
    temp_table_multi_final = pd.concat([temp_table_multi_clean, temp_table_multi_chimeric], ignore_index=True)
    Def_barcodes = temp_table_multi_final[temp_table_multi_final['Mode'] == 'Def']
    Chimeric_barcodes = temp_table_multi_final[temp_table_multi_final['Mode'] == 'Chimeric']
    
    # Get the number of unique barcodes in each table
    num_unique_clean = len(temp_table_multi_clean['BC'].unique())
    num_unique_chimeric = len(temp_table_multi_chimeric['BC'].unique())
    num_unique_single = len(temp_table_single['BC'].unique())
    total_barcodes = num_unique_clean + num_unique_chimeric + num_unique_single
    
    logger.info(f"Number of single read barcodes mapping to only one fragment: Reads:{len(temp_table_single)} Barcodes:{num_unique_single} ({num_unique_single/total_barcodes*100:.2f}%)")
    logger.info(f"Number of multi read barcodes mapping to only one fragment: Reads:{len(temp_table_multi_clean)} Barcodes:{num_unique_clean} ({num_unique_clean/total_barcodes*100:.2f}%)")
    logger.info(f"   Of theses are cleaned chimeric barcodes (ratio above {threshold}) Reads:{len(temp_table_multi_chimeric[temp_table_multi_chimeric['Mode'] == 'Def'])} Barcodes: {len(temp_table_multi_chimeric[temp_table_multi_chimeric['Mode'] == 'Def']['BC'].unique())} ({len(temp_table_multi_chimeric[temp_table_multi_chimeric['Mode'] == 'Def']['BC'].unique())/total_barcodes*100:.2f}%)")
    logger.info(f"Number of multi read barcodes mapping to multiple fragments: Reads:{len(temp_table_multi_chimeric)} Barcodes:{num_unique_chimeric} ({num_unique_chimeric/total_barcodes*100:.2f}%)")

    Def_barcodes = temp_table_multi_final[temp_table_multi_final['Mode'] == 'Def']
    Chimeric_barcodes = temp_table_multi_final[temp_table_multi_final['Mode'] == 'Chimeric']
    
    return Def_barcodes, Chimeric_barcodes, temp_table_single


def write_def_barcodes(Def_barcodes: pd.DataFrame, out_name: str, name_suffix:str="")-> None:
    """
    Write the definitiv barcodes to a CSV file.
    
    Parameters:
        Def_barcodes (pd.DataFrame): DataFrame containing the definitiv barcodes
        out_name (str): Path to the output file
        
    Returns:
        None
    """
    # split out_name to get the path
    out_path = "/".join(out_name.split("/")[:-1])
    save_path = out_path + f"/barcode_db{name_suffix}.fasta"
    # remove the file if it already exists
    if os.path.exists(save_path):
        os.remove(save_path)
    # get all unique barcodes from the Def_barcodes table
    Def_barcodes = Def_barcodes.drop_duplicates(subset=['BC'])
    # save the definitiv barcodes from the BC column to a fasta file
    for i, row in Def_barcodes.iterrows():
        with open(save_path, 'a') as f:
            f.write(f">BC_{i+1}\n{row['BC']}\n")


def main():
    start_time = datetime.now()

    # load config
    config = get_config("S3")
    
    # Create a logger
    create_logger(config["log_dir"], "S3")

    # read in the LUT file
    lut_df = pd.read_csv(config["in_name_LUT"])
    
    # Save unique fragments and barcodes as FASTA file
    save_unique_fragments_barcodes(config["fragment_file"], config["barcode_file"])
    
    # output file_name
    output_file = config["out_name"].replace(".csv", ".h5")
    
    chunk_size = config["chunk_size"]
    write_mode = 'w'
    
    # Process fragments and barcodes in chunks
    for i, (frag_chunk, bc_chunk) in enumerate(load_frag_bc_reads_chunked(config["fragment_file"], config["barcode_file"], chunk_size)):

        # Create full table for the chunk
        chunk_table = create_full_table(lut_df, frag_chunk, bc_chunk)

        # Write chunk directly to final output file in HDF5 or Parquet format
        chunk_table.to_hdf(output_file, key='data', mode=write_mode, format='table', append=(write_mode == 'a'))

        # After the first write, change mode to 'append'
        write_mode = 'a'

        # Explicitly free memory
        del frag_chunk, bc_chunk, chunk_table
        gc.collect()
    
    # read the full table from the output file
    full_table = pd.read_hdf(output_file, key='data')
    
    # remove the output file
    os.remove(output_file)

    # Perform Starcode barcode reduction and replace barcodes with reduced versions
    full_table = starcode_based_reduction_and_replace(full_table, config['barcode_file'], 'BC')

    # Split reads into single-read and multi-read barcodes
    temp_table_single, temp_table_multi = split_reads_into_single_and_multi_read_barcodes(full_table)
    #remove full_table df from memory
    del full_table

    # Split multi-read barcodes into clean and chimeric
    temp_table_multi_clean, temp_table_multi_chimeric = split_multi_read_barcodes_into_clean_and_chimeric(temp_table_multi)

    # Create consensus alignment for chimeric barcodes and get valid chimeric barcodes based on threshold
    temp_table_multi_chimeric = get_valid_chimeric_barcodes(temp_table_multi_chimeric, config["threshold"])

    # Combine all tables into final output
    def_barcodes_table, chimeric_barcode_table, single_barcode_table = combine_tables(temp_table_multi_clean, temp_table_multi_chimeric, temp_table_single, config["threshold"])
    del temp_table_multi_clean, temp_table_multi_chimeric, temp_table_single
    
    # Add the single barcodes to the final table or chimeric barcodes
    final_barcodes_table = def_barcodes_table.copy()
    if config["single_read"]:
        # merge single read barcodes with the multi read barcodes
        final_barcodes_table = pd.concat([final_barcodes_table, single_barcode_table], ignore_index=True)
    if config["chimeric_read"]:
        # merge chimeric barcodes with the multi read barcodes
        final_barcodes_table = pd.concat([final_barcodes_table, chimeric_barcode_table], ignore_index=True)

    # Save the output tables
    final_barcodes_table.to_csv(config['out_name'], index=False)
    logger.info(f"Final barcodes saved to {config['out_name']}")
    def_barcodes_table.to_csv(config['out_name'].replace(".csv", "_def.csv"), index=False)
    logger.info(f"Definitiv barcodes saved to {config['out_name']}")
    chimeric_barcode_table.to_csv(config['out_name'].replace(".csv", "_chimeric.csv"), index=False)
    logger.info(f"Chimeric barcodes saved to {config['out_name'].replace('.csv', '_chimeric.csv')}")
    single_barcode_table.to_csv(config['out_name'].replace(".csv", "_single.csv"), index=False)
    logger.info(f"Single barcodes saved to {config['out_name'].replace('.csv', '_single.csv')}")
    
    allowed = "Definitiv"
    if config["single_read"]:
        allowed += ", Single"
    if config["chimeric_read"]:
        allowed += ", Chimeric"
        
    logger.info(f"Included barcodes: {allowed}")
    
    # write the definitiv barcodes to a fasta file
    write_def_barcodes(final_barcodes_table, config['out_name'])
    
    # Print total analysis time
    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
