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
        BC,LUTnr,bitScore,mismatches,tCount,mCount,Mode
"""



import psutil
import gzip
import os
import tempfile
import subprocess
import multiprocessing
import pandas as pd
from Bio import SeqIO
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


def run_command(command: list, description: str) -> tuple:
    """
    Runs a subprocess command and returns stdout, stderr, and error status.
    
    Parameters: 
        command (list): The command to run
        description (str): A description of the command
    
    Returns:   
        tuple: stdout and stderr
    """
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            logger.error(f"Error running {description} with code {process.returncode}")
            logger.error(stderr)
            sys.exit(1)

        return stdout, stderr


def load_frag_bc_reads_chunked(fragments_file: str, barcodes_file: str, memory_fraction: float = 0.8):
    """
    Load the fragments and barcodes from FASTQ files in dynamically sized chunks based on memory.

    Parameters:
        fragments_file (str): Path to the fragments FASTQ file.
        barcodes_file (str): Path to the barcodes FASTQ file.
        memory_fraction (float): Fraction of available memory to use (default: 0.8).

    Yields:
        tuple: Lists of fragments and barcodes as SeqRecords.
    """
    def get_optimal_chunk_size(average_record_size: int):
        # Get total available memory and allocate per thread
        available_memory = psutil.virtual_memory().available * memory_fraction
        num_threads = multiprocessing.cpu_count()
        memory_per_thread = available_memory / num_threads
        
        # Calculate chunk size based on memory per thread and average record size
        return int(memory_per_thread // average_record_size)

    def estimate_average_record_size(file_path):
        # Sample a few records to estimate their average size in memory
        with gzip.open(file_path, "rt") as handle:
            records = list(SeqIO.parse(handle, "fastq"))
            if not records:
                return 0
            return sum(len(record.seq) + len(record.letter_annotations.get("phred_quality", []))
                       for record in records) // len(records)

    # Estimate average record sizes for fragments and barcodes
    avg_frag_size = estimate_average_record_size(fragments_file)
    avg_bc_size = estimate_average_record_size(barcodes_file)

    # Use the larger average record size for a conservative estimate
    avg_record_size = max(avg_frag_size, avg_bc_size)

    # Get the optimal chunk size
    chunk_size = get_optimal_chunk_size(avg_record_size)

    def read_fastq_chunks(file_path):
        with gzip.open(file_path, "rt") as handle:
            chunk = []
            for record in SeqIO.parse(handle, "fastq"):
                chunk.append(record)
                if len(chunk) >= chunk_size:
                    yield chunk
                    chunk = []
            if chunk:
                yield chunk

    frag_gen = read_fastq_chunks(fragments_file)
    bc_gen = read_fastq_chunks(barcodes_file)

    for frag_chunk, bc_chunk in zip(frag_gen, bc_gen):
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

    return full_table


def create_LUTnrs(full_table: pd.DataFrame) -> pd.DataFrame:
    """
    Create a DataFrame with LUTnrs and their corresponding fragments.

    Parameters:
        full_table (pd.DataFrame): DataFrame containing the full table with fragments and barcodes
        
    Returns:
        pd.DataFrame: DataFrame containing the LUTnrs a for every fragment
    """
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

    return lutnrs_df


def starcode_based_barcode_reduction(full_table: pd.DataFrame) -> pd.DataFrame:
    """
    Perform barcode reduction using Starcode clustering on the unique barcodes from the full table.

    Parameters:
        full_table (pd.DataFrame): DataFrame containing the full table with barcodes (column 'BC')
        
    Returns:
        pd.DataFrame: DataFrame containing the Starcode-reduced barcodes
    """
    logger.info("Running Starcode clustering on unique barcodes from the full table")
    
    # Extract unique barcodes
    unique_barcodes = full_table['BC'].unique()
    
    # Write unique barcodes to a temporary file
    barcode_temp_file = tempfile.NamedTemporaryFile(delete=False, mode='w+')
    barcode_temp_file.writelines(f"{bc}\n" for bc in unique_barcodes)
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

    logger.info("Starcode clustering completed")
    return starcode_exploded


def replace_barcodes_with_starcode_versions(full_table: pd.DataFrame, starcode_exploded: pd.DataFrame)-> pd.DataFrame:
    """
    Replace barcodes in the full table with their Starcode-reduced versions.

    Parameters:
        full_table (pd.DataFrame): DataFrame containing the full table with fragments and barcodes
        starcode_exploded (pd.DataFrame): DataFrame containing the Starcode-reduced barcodes
        
    Returns:
        pd.DataFrame: DataFrame containing the full table with Starcode-reduced barcodes
    """
    logger.info("Replacing barcodes with Starcode-reduced versions")
    # combine the full table with the starcode exploded table based on the BC
    full_table = full_table.merge(starcode_exploded[['BC', 'scBC']], on='BC', how='left')
    # rename the coulmns and drop the old BC column
    full_table.rename(columns={'BC': 'oldBC', 'scBC': 'BC'}, inplace=True)
    full_table.drop(columns=['oldBC'], inplace=True)

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
    # if their is only one BC for a given LUTnr, then it is a clean barcode
    clean_barcodes = bc_counts[bc_counts == 1].index
    # if their is more than one BC for a given LUTnr, then it is a chimeric barcode
    chimeric_barcodes = bc_counts[bc_counts > 1].index

    # create tables with clean and chimeric barcodes
    temp_table_multi_clean = temp_table_multi_grouped[temp_table_multi_grouped['BC'].isin(clean_barcodes)].copy()
    temp_table_multi_chimeric = temp_table_multi_grouped[temp_table_multi_grouped['BC'].isin(chimeric_barcodes)].copy()

    # set the max read count, total read count, since we only have one barcode for this LUTnr and mode for clean barcode
    temp_table_multi_clean['mCount'] = temp_table_multi_clean['tCount']
    temp_table_multi_clean['Mode'] = 'Def'

    return temp_table_multi_clean, temp_table_multi_chimeric


def calculate_consensus_alignment(temp_table_multi_chimeric: pd.DataFrame)-> pd.DataFrame:
    """
    Calculate consensus alignment of chimeric barcodes. This means we are getting the barcode with the highest maximal read count for each LUTnr.

    Parameters:
        temp_table_multi_chimeric (pd.DataFrame): DataFrame containing the chimeric multi-read barcodes
    
    Returns:
        pd.DataFrame: DataFrame containing the consensus alignment for chimeric barcodes
    """
    logger.info("Calculating consensus alignment for chimeric barcodes")
    # For every barcode LUTnr pair, set mCount to tCount
    temp_table_multi_chimeric['mCount'] = temp_table_multi_chimeric['tCount']
    # Sum the total count of reads for each barcode
    temp_table_multi_chimeric['tCount'] = temp_table_multi_chimeric.groupby('BC')['tCount'].transform('sum')

    # Select the LUTnr, barcode pair with the highest maximal read count
    idx = temp_table_multi_chimeric.groupby('BC')['mCount'].idxmax()
    temp_table_multi_consensus = temp_table_multi_chimeric.loc[idx].copy()
    temp_table_multi_consensus['Mode'] = 'Def'

    return temp_table_multi_consensus


def combine_tables(temp_table_multi_clean: pd.DataFrame, temp_table_multi_consensus: pd.DataFrame, temp_table_single: pd.DataFrame)-> pd.DataFrame:
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
    logger.info(f"Number of barcodes mapping to only one fragment: {len(temp_table_multi_clean)}")
    logger.info(f"Number of barcodes mapping to more than one fragment: {len(temp_table_multi_consensus)}")
    logger.info(f"Number of barcodes-fragment pairs sequenced only once: {len(temp_table_single)}")
    output_table = pd.concat([temp_table_multi_final, temp_table_single], ignore_index=True)
    # remove the reads column
    output_table.drop(columns=['Reads'], inplace=True)

    return output_table

def process_chunk(frag_chunk, bc_chunk):
        """
        Process a single chunk of fragments and barcodes.

        Parameters:
            frag_chunk (list): Chunk of fragment reads as SeqRecords
            bc_chunk (list): Chunk of barcode reads as SeqRecords
            lut_df (pd.DataFrame): LUT DataFrame
            blast_db_prefix (str): BLAST database prefix

        Returns:
            pd.DataFrame: Full table for the processed chunk
        """

        # Create full table with fragments and barcodes
        chunk_table = create_full_table(frag_chunk, bc_chunk)
 
        return chunk_table


def main():
    start_time = datetime.now()

    # load config
    config = get_config("S3")
    
    # Create a logger
    create_logger(config["log_dir"], "S3")

    # Parallelize processing of chunks
    chunk_results = []
    chunk_processing_args = []

    # Prepare arguments for each chunk
    for frag_chunk, bc_chunk in load_frag_bc_reads_chunked(config["fragment_file"], config["barcode_file"]):
        chunk_processing_args.append((frag_chunk, bc_chunk))
    
    logger.info(f"Number of chunks to process: {len(chunk_processing_args)}")

    # Use multiprocessing Pool to process chunks in parallel
    with multiprocessing.Pool(processes=multiprocessing.cpu_count() - 1) as pool:
        chunk_results = pool.starmap(process_chunk, chunk_processing_args)
    
    # Concatenate all chunk tables
    full_table = pd.concat(chunk_results, ignore_index=True)
    # log the number of unique fragments and barcodes
    logger.info(f"Number of unique fragments: {len(full_table['Reads'].unique())}")
    logger.info(f"Number of unique barcodes: {len(full_table['BC'].unique())}")

    # Create LUTnrs
    LUT_df = create_LUTnrs(full_table)
        
    # Merge the LUTnrs with the full table but only add the LUTnr column
    full_table = full_table.merge(LUT_df[["Reads","LUTnr"]], on='Reads', how='left')
    
    # Perform Starcode barcode reduction
    starcode_exploded = starcode_based_barcode_reduction(full_table)

    # Replace barcodes with Starcode-reduced versions
    full_table = replace_barcodes_with_starcode_versions(full_table, starcode_exploded)

    # Split reads into single-read and multi-read barcodes
    temp_table_single, temp_table_multi = split_reads_into_single_and_multi_read_barcodes(full_table)

    # Split multi-read barcodes into clean and chimeric
    temp_table_multi_clean, temp_table_multi_chimeric = split_multi_read_barcodes_into_clean_and_chimeric(temp_table_multi)

    # Calculate consensus alignment for chimeric barcodes
    temp_table_multi_consensus = calculate_consensus_alignment(temp_table_multi_chimeric)

    # Combine all tables into final output
    output_table = combine_tables(temp_table_multi_clean, temp_table_multi_consensus, temp_table_single)
    
    # add the LUTnr peptide to the output table
    output_table = output_table.merge(LUT_df, on="LUTnr", how='left')
    
    # change reads to Sequence
    output_table.rename(columns={'Reads': 'Sequence'}, inplace=True)

    # Save the output table
    output_table.to_csv(config['out_name'], index=False)
    logger.info(f"Output saved to {config['out_name']}")

    # Print total analysis time
    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
