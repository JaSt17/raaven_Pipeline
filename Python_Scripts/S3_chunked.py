#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script maps the found fragments and barcodes to the LUT data.
It then reduces the barcodes using the Starcode algorithm. And replaces the barcodes with the Starcode-reduced versions.
Finally it devides the barcodes into single-read and multi-read barcodes and splits the multi-read barcodes into clean and chimeric barcodes.
Saves all found barcodes to a CSV file.

Workflow:
    - Load the LUT data
    - Create a BLAST database from the LUT sequences
    - Load fragments and barcodes from FASTQ files in chunks
        - For each chunk:
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

import gzip
import os
import tempfile
import subprocess
import multiprocessing
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from datetime import datetime
from itertools import islice
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


def load_frag_bc_reads_chunked(fragments_file: str, barcodes_file: str, num_chunks: int):
    """
    Load the fragments and barcodes from FASTQ files, splitting the data into equally sized parts based on the number of chunks.

    Parameters:
        fragments_file (str): Path to the fragments FASTQ file.
        barcodes_file (str): Path to the barcodes FASTQ file.
        num_chunks (int): Number of chunks to split the data into.

    Yields:
        tuple: Lists of fragments and barcodes as SeqRecords.
    """

    def count_fastq_records(file_path):
        """
        Counts the number of records in a gzipped FASTQ file using zcat and wc -l.

        Parameters:
            file_path (str): Path to the FASTQ file.

        Returns:
            int: Total number of records.
        """
        cmd = [f"zcat {file_path} | wc -l"]
        description = "Counting lines in FASTQ file"

        stdout, stderr = run_command(cmd, description, shell=True)

        if stderr:
            logger.error(f"Error counting lines: {stderr}")
            raise Exception(f"Error counting lines: {stderr}")

        line_count = int(stdout.strip())
        total_records = line_count // 4
        return total_records

    total_frag_records = count_fastq_records(fragments_file)
    total_bc_records = count_fastq_records(barcodes_file)

    if total_frag_records != total_bc_records:
        logger.error("Number of records in fragments and barcodes files do not match.")
        raise ValueError("Number of records in fragments and barcodes files do not match.")

    total_records = total_frag_records
    logger.info(f"Total number of barcode-fragment-pairs: {total_records}")

    # Calculate chunk sizes
    base_chunk_size = total_records // num_chunks
    remainder = total_records % num_chunks

    # Create a list of chunk sizes that distribute the remainder evenly
    chunk_sizes = [base_chunk_size + 1 if i < remainder else base_chunk_size for i in range(num_chunks)]

    def read_fastq_chunks(file_path):
        """
        Generator function to read chunks from a FASTQ file.

        Parameters:
            file_path (str): Path to the FASTQ file.

        Yields:
            list: A chunk of SeqRecords.
        """
        with gzip.open(file_path, "rt") as handle:
            record_iter = SeqIO.parse(handle, "fastq")
            for chunk_size in chunk_sizes:
                chunk = list(islice(record_iter, chunk_size))
                if chunk:
                    yield chunk
                else:
                    break

    frag_gen = read_fastq_chunks(fragments_file)
    bc_gen = read_fastq_chunks(barcodes_file)

    for frag_chunk, bc_chunk in zip(frag_gen, bc_gen):
        yield frag_chunk, bc_chunk


def make_customarray_reference_index(lut_df: pd.DataFrame)-> str:
    """
    Create a BLAST database from LUT sequences.

    Parameters:
        lut_df (pd.DataFrame): DataFrame containing the LUT sequences

    Returns:
        str: Path to the BLAST database prefix
    """
    logger.info("Creating BLAST database from LUT sequences")
    lut_fa = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa")
    lut_records = [SeqRecord(Seq(row['Sequence']), id=str(row['LUTnr']), description="")
                    for _, row in lut_df.iterrows()]
    SeqIO.write(lut_records, lut_fa.name, "fasta")
    lut_fa.close()

    # Create BLAST database
    blast_db_prefix = tempfile.mktemp(prefix="blastDB_")
    makeblastdb_cmd = f"makeblastdb -in {lut_fa.name} -out {blast_db_prefix} -dbtype nucl -title LUT -parse_seqids -logfile /dev/null"
    stdout, stderr = run_command(makeblastdb_cmd.split(), "BLAST database creation")
    # print the size of the database
    check_db_size_cmd = f"blastdbcmd -db {blast_db_prefix} -info"
    stdout, stderr = run_command(check_db_size_cmd.split(), "Checking BLAST database size")
    # pritn the first 2 lines of the output
    stdout = stdout.split('\n')[0:2]
    size = stdout[1].replace("\t", "")
    logger.info(f"BLAST database size: {size}")
    return blast_db_prefix


def save_unique_fragments(reads_frag: list)-> tuple:
    """
    Save unique all found unique fragments from the library to a FASTA file.

    Parameters:
        reads_frag (list): List of fragments as SeqRecords
        
    Returns:
        tuple: Path to the FASTA file containing unique fragments and the unique fragments
    """
    unique_fragments = set(str(rec.seq) for rec in reads_frag)
    
    # exclude all sequences that have a N in them
    unique_fragments = [seq for seq in unique_fragments if 'N' not in seq]
    
    unique_reads = [SeqRecord(Seq(seq), id=str(i+1), description="")
                    for i, seq in enumerate(unique_fragments)]
    fragments_unique_fa = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa")
    SeqIO.write(unique_reads, fragments_unique_fa.name, "fasta-2line")
    fragments_unique_fa.close()
    return fragments_unique_fa.name, unique_fragments


def align_against_library(fragments_unique_fa_name: str, blast_db_prefix: str)-> str:
    """
    Align unique fragments against the LUT database using BLASTn. BLAST output is saved to a file.

    Parameters:
        fragments_unique_fa_name (str): Path to the FASTA file containing unique fragments
        blast_db_prefix (str): Path to the BLAST database prefix
    
    Returns:
        str: Path to the BLAST output file
    """
    # Get the number of threads
    num_threads = multiprocessing.cpu_count()
    blast_out_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".txt")
    blast_cmd = [
        "blastn",
        "-query", fragments_unique_fa_name,
        "-db", blast_db_prefix,
        "-outfmt", "10",
        "-max_target_seqs", "25",
        "-word_size", "7",
        "-num_threads", str(num_threads),
        "-out", blast_out_file.name
    ]
    stdout, stderr =  run_command(blast_cmd, "BLAST alignment")

    # Check for errors
    if stderr:
        logger.error(f"BLAST error: {stderr}")
        raise Exception(f"BLAST error: {stderr}")
    
    return blast_out_file.name


def read_blast_output(blast_out_file_name: str, unique_fragments: set, lut_df: pd.DataFrame)-> pd.DataFrame:
    """
    Read the BLAST output into a DataFrame and map every read to its corresponding LUTnr.
    
    Parameters:
        blast_out_file_name (str): Path to the BLAST output file
        unique_fragments (set): Set of unique fragments
        lut_df (pd.DataFrame): DataFrame containing the LUT sequences
    
    Returns:
        pd.DataFrame: DataFrame containing the BLAST results
    """
    blast_columns = ["Reads", "Sequence", "identity", "alignmentLength", "mismatches",
                    "gapOpens", "q_start", "q_end", "s_start", "s_end",
                    "evalue", "bitScore"]
    blast_df = pd.read_csv(blast_out_file_name, header=None, names=blast_columns)
    
    # Map the Reads and Sequence IDs back to the actual sequences
    reads_index = {str(i+1): seq for i, seq in enumerate(unique_fragments)}
    lut_seq_index = lut_df.set_index('LUTnr')['Sequence'].astype(str).to_dict()
    blast_df['Reads'] = blast_df['Reads'].astype(str).map(reads_index)
    blast_df['Sequence'] = blast_df['Sequence'].astype(str).map(lut_seq_index)
    
    return blast_df


def create_full_table(blast_df: pd.DataFrame, lut_df: pd.DataFrame, reads_frag: list, reads_BC: list)-> pd.DataFrame:
    """
    Create a full table of all fragments that matched the LUT with their BLASTn results. For every fragment, only the top hit is kept.
    
    Parameters:
        blast_df (pd.DataFrame): DataFrame containing the BLAST results
        lut_df (pd.DataFrame): DataFrame containing the LUT sequences
        reads_frag (list): List of fragments as SeqRecords
        reads_BC (list): List of barcodes as SeqRecords
    
    Returns:
        pd.DataFrame: DataFrame containing the full table with fragments and barcodes
    """
    # Merge BLAST results with LUT data
    blast_df = blast_df.merge(lut_df[['Sequence', 'LUTnr']], on='Sequence', how='inner')
    blast_df.drop(columns=['Sequence'], inplace=True)

    # Convert columns to appropriate data types
    blast_df['bitScore'] = pd.to_numeric(blast_df['bitScore'])
    blast_df['mismatches'] = pd.to_numeric(blast_df['mismatches'])

    # Sort by Reads, LUTnr, and bitScore (descending for bitScore)
    blast_df.sort_values(by=['Reads', 'LUTnr', 'bitScore'], ascending=[True, True, False], inplace=True)
    # Keep only the first occurrence of each unique pair of Reads and LUTnr
    blast_df.drop_duplicates(subset=['Reads', 'LUTnr'], keep='first', inplace=True)

    # Create full table with fragments and barcodes pairs
    full_table = pd.DataFrame({
        'Reads': [str(rec.seq) for rec in reads_frag],
        'BC': [str(rec.seq) for rec in reads_BC]
    })

    # Selecting only the top hit for each fragment from the BLAST results
    blast_top_hit = blast_df.loc[blast_df.groupby('Reads')['bitScore'].idxmax()]
    full_table = full_table.merge(blast_top_hit, on='Reads', how='inner')
    
    # Calculate the percentage of fragments that have been aligned to our DB (original created fragments)
    alignment_percentage = len(full_table) / len(reads_frag) if reads_frag else 0
    logger.info(f"Percentage of fragments aligned to the LUT: {alignment_percentage:.2%}")

    return full_table


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
    temp_table_multi['mismatches'] = temp_table_multi['mismatches'].astype(float)

    # Group by BC and LUTnr to compute statistics for every single pair
    temp_table_multi_grouped = temp_table_multi.groupby(['BC', 'LUTnr']).agg({
        'bitScore': 'mean',
        'mismatches': 'median',
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

    return output_table

def process_chunk(args):
    """
    Process a single chunk of fragments and barcodes.

    Parameters:
        args (tuple): Tuple containing the following elements:
            frag_chunk (list): Chunk of fragment reads as SeqRecords
            bc_chunk (list): Chunk of barcode reads as SeqRecords
            lut_df (pd.DataFrame): LUT DataFrame
            blast_db_prefix (str): BLAST database prefix
    """
    frag_chunk, bc_chunk, lut_df, blast_db_prefix = args

    # Proceed with processing as before
    frag_fa_name, unique_fragments = save_unique_fragments(frag_chunk)
    blast_out_file_name = align_against_library(frag_fa_name, blast_db_prefix)
    blast_df = read_blast_output(blast_out_file_name, unique_fragments, lut_df)
    chunk_table = create_full_table(blast_df, lut_df, frag_chunk, bc_chunk)
    
    # Clean up temporary fragment FASTA file
    try:
        os.remove(frag_fa_name)
        os.remove(blast_out_file_name)
    except Exception as e:
        logger.warning(f"Failed to delete temporary files: {e}")
    
    return chunk_table


def main():
    start_time = datetime.now()

    # Load config
    config = get_config("S3")

    # Create main logger
    create_logger(config["log_dir"], "S3")

    # Read LUT file
    lut_df = pd.read_csv(config["in_name_LUT"])

    # Create BLAST database
    blast_db_prefix = make_customarray_reference_index(lut_df)
    
    # Determine the number of worker processes
    num_workers = 5
    
    # Prepare arguments for each chunk
    chunk_args = []
    for frag_chunk, bc_chunk in load_frag_bc_reads_chunked(config["fragment_file"], config["barcode_file"], num_workers):
        chunk_args.append((frag_chunk, bc_chunk, lut_df, blast_db_prefix))
    
    #create a temporary directory to store the chunk files
    temp_dir = tempfile.mkdtemp()
    
    # Initialize multiprocessing Pool
    with multiprocessing.Pool(processes=min(10, num_workers)) as pool:
        # Use imap_unordered for better performance
        chunk_results = pool.imap_unordered(process_chunk, chunk_args)
        for i, chunk_table in enumerate(chunk_results):
            chunk_table.to_pickle(temp_dir+f"/chunk_{i}.pkl")
    
    # Load all chunk tables
    chunk_tables = [pd.read_pickle(temp_dir+f"/chunk_{i}.pkl") for i in range(num_workers)]
    
    # Concatenate all chunk tables
    full_table = pd.concat(chunk_tables, ignore_index=True)

    # Log the number of unique fragments and barcodes
    logger.info(f"Number of unique fragments: {full_table['Reads'].nunique()}")
    logger.info(f"Number of unique barcodes: {full_table['BC'].nunique()}")

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

    # Save the output table
    output_table.to_csv(config['out_name'], index=False)
    logger.info(f"Output saved to {config['out_name']}")
    logger.info(f"Total execution time: {datetime.now() - start_time}")


if __name__ == "__main__":
    main()
