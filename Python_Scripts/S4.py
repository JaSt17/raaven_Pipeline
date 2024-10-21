#!/usr/bin/env python3
"""
Author: Jaro Steindorff
"""

import os
import gzip
import tempfile
import subprocess
import multiprocessing
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import time
import logging
# local import
from config import get_config

# Initialize logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def load_trimmed_reads(fragments_file: str, barcodes_file: str)-> tuple:
    """
    Load the trimmed reads and barcodes from FASTQ files.

    :param fragments_file: Path to the trimmed reads FASTQ file
    :param barcodes_file: Path to the barcodes FASTQ file

    :return: A tuple containing the trimmed reads and barcodes as lists of SeqRecords
    """
    logger.info("Reading trimmed reads and barcodes")
    with gzip.open(fragments_file, "rt") as handle:
        reads_trim = list(SeqIO.parse(handle, "fastq"))
    with gzip.open(barcodes_file, "rt") as handle:
        reads_BC = list(SeqIO.parse(handle, "fastq"))
    return reads_trim, reads_BC


def make_customarray_reference_index(lut_df: pd.DataFrame)-> str:
    """
    Create a BLAST database from LUT sequences.

    :param lut_df: DataFrame containing the LUT sequences

    :return: Path to the BLAST database prefix
    """
    logger.info("Creating BLAST database from LUT sequences")
    lut_fa = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa")
    lut_records = [SeqRecord(Seq(row['Sequence']), id=str(row['LUTnr']), description="")
                    for _, row in lut_df.iterrows()]
    SeqIO.write(lut_records, lut_fa.name, "fasta")
    lut_fa.close()

    # Create BLAST database
    blast_db_prefix = tempfile.mktemp(prefix="blastDB_")
    makeblastdb_cmd = f"makeblastdb -in {lut_fa.name} -out {blast_db_prefix} -dbtype nucl -title LUT -parse_seqids"
    subprocess.run(makeblastdb_cmd, shell=True, check=True)
    return blast_db_prefix


def save_unique_fragments(reads_trim: list)-> tuple:
    """
    Save unique fragments from trimmed reads to a FASTA file.

    :param reads_trim: List of trimmed reads as SeqRecords

    :return: Path to the FASTA file containing unique fragments and a set of unique sequences
    """
    logger.info("Saving unique fragments to FASTA file")
    unique_reads_sequences = set(str(rec.seq) for rec in reads_trim)
    unique_reads = [SeqRecord(Seq(seq), id=str(i+1), description="")
                    for i, seq in enumerate(unique_reads_sequences)]
    fragments_unique_fa = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".fa")
    SeqIO.write(unique_reads, fragments_unique_fa.name, "fasta")
    fragments_unique_fa.close()
    return fragments_unique_fa.name, unique_reads_sequences


def align_against_library(fragments_unique_fa_name: str, blast_db_prefix: str)-> str:
    """
    Align unique fragments against the LUT database using BLASTn.

    :param fragments_unique_fa_name: Path to the FASTA file containing unique fragments
    :param blast_db_prefix: Path to the BLAST database prefix

    :return: Path to the BLAST output file
    """
    logger.info("Running BLAST alignment")
    # Get the number of threads
    num_threads = multiprocessing.cpu_count()
    blast_out_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".txt")
    blastn_cline = NcbiblastnCommandline(query=fragments_unique_fa_name,
                                        db=blast_db_prefix,
                                        outfmt=10,
                                        max_target_seqs=25,
                                        word_size=7,
                                        num_threads=num_threads,
                                        out=blast_out_file.name)
    stdout, stderr = blastn_cline()

    # Check for errors
    if stderr:
        logger.error(f"BLAST error: {stderr}")
        raise Exception(f"BLAST error: {stderr}")
    
    return blast_out_file.name


def compress_blast_output(blast_out_file_name: str, output_filename: str)-> None:
    """
    Compress the BLAST output and save it.

    :param blast_out_file_name: Path to the BLAST output file
    :param output_filename: Path to the compressed output file

    :return: None
    """
    logger.info("Compressing BLAST output")
    with open(blast_out_file_name, 'rb') as f_in, gzip.open(output_filename, 'wb') as f_out:
        f_out.writelines(f_in)


def read_blast_output(blast_out_file_name: str, unique_reads_sequences: set, lut_df: pd.DataFrame)-> pd.DataFrame:
    """
    Read the BLAST output into a DataFrame and map IDs to sequences.
    """
    logger.info("Reading BLAST output")
    blast_columns = ["Reads", "Sequence", "identity", "alignmentLength", "mismatches",
                    "gapOpens", "q_start", "q_end", "s_start", "s_end",
                    "evalue", "bitScore"]
    blast_df = pd.read_csv(blast_out_file_name, header=None, names=blast_columns)
    
    # Map the Reads and Sequence IDs back to the actual sequences
    reads_index = {str(i+1): seq for i, seq in enumerate(unique_reads_sequences)}
    lut_seq_index = lut_df.set_index('LUTnr')['Sequence'].astype(str).to_dict()
    blast_df['Reads'] = blast_df['Reads'].astype(str).map(reads_index)
    blast_df['Sequence'] = blast_df['Sequence'].astype(str).map(lut_seq_index)
    
    return blast_df


def create_full_table(blast_df: pd.DataFrame, lut_df: pd.DataFrame, reads_trim: list, reads_BC: list)-> pd.DataFrame:
    """
    Create a full table of all fragments that matched the LUT with their BLASTn results.
    
    :param blast_df: DataFrame containing the BLAST results
    :param lut_df: DataFrame containing the LUT sequences
    :param reads_trim: List of trimmed reads as SeqRecords
    :param reads_BC: List of barcodes as SeqRecords
    """
    logger.info("Creating full table with fragments and barcodes")
    # Merge BLAST results with LUT data
    blast_df = blast_df.merge(lut_df[['Sequence', 'LUTnr']], on='Sequence', how='left')
    blast_df.dropna(subset=['LUTnr'], inplace=True)
    blast_df.drop(columns=['Sequence'], inplace=True)

    # Convert columns to appropriate data types
    blast_df['bitScore'] = pd.to_numeric(blast_df['bitScore'])
    blast_df['mismatches'] = pd.to_numeric(blast_df['mismatches'])

    # Keep only the top hit for each read and LUTnr
    blast_df.sort_values(by=['Reads', 'LUTnr', 'bitScore'], ascending=[True, True, False], inplace=True)
    blast_df.drop_duplicates(subset=['Reads', 'LUTnr'], keep='first', inplace=True)

    # Create full table with fragments and barcodes
    full_table = pd.DataFrame({
        'Reads': [str(rec.seq) for rec in reads_trim],
        'BC': [str(rec.seq) for rec in reads_BC]
    })

    # Selecting only the top hit for each read
    blast_top_hit = blast_df.loc[blast_df.groupby('Reads')['bitScore'].idxmax()]
    full_table = full_table.merge(blast_top_hit, on='Reads', how='inner')

    # Calculate the percentage of reads that have been aligned to our DB
    alignment_percentage = len(full_table) / len(reads_trim) if reads_trim else 0
    logger.info(f"Alignment percentage: {alignment_percentage:.2%}")

    return full_table


def starcode_based_barcode_reduction(barcodes_file: str)-> pd.DataFrame:
    """
    Perform barcode reduction using Starcode clustering.

    :param barcodes_file: Path to the barcode file

    :return: DataFrame containing the Starcode-reduced barcodes
    """
    logger.info("Running Starcode clustering")
    # get the number of threads
    num_threads = multiprocessing.cpu_count()
    starcode_output = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix=".txt")
    starcode_cmd = f"gzip -cd {barcodes_file} | starcode -t {num_threads-1} --print-clusters -d 1 -r5 -q -o {starcode_output.name}"
    subprocess.run(starcode_cmd, shell=True, check=True)

    # Read Starcode output
    starcode_df = pd.read_csv(starcode_output.name, sep='\t', header=None, names=['scBC', 'count', 'BC_list'])
    starcode_df['BC_list'] = starcode_df['BC_list'].str.split(',')

    # Explode the BC_list to have one row per barcode
    starcode_exploded = starcode_df.explode('BC_list')
    starcode_exploded.rename(columns={'BC_list': 'BC'}, inplace=True)

    return starcode_exploded


def replace_barcodes_with_starcode_versions(full_table: pd.DataFrame, starcode_exploded: pd.DataFrame)-> pd.DataFrame:
    """
    Replace barcodes in the full table with their Starcode-reduced versions.

    :param full_table: DataFrame containing the full table with fragments and barcodes
    :param starcode_exploded: DataFrame containing the Starcode-reduced barcodes

    :return: DataFrame containing the full table with Starcode-reduced barcodes
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

    :param full_table: DataFrame containing the full table with fragments and barcodes

    :return: Tuple containing the single-read and multi-read tables
    """
    logger.info("Splitting reads into single-read and multi-read barcodes")
    barcode_counts = full_table['BC'].value_counts()
    # get single read barcodes table
    single_read_barcodes = barcode_counts[barcode_counts == 1].index
    # get multi read barcodes table
    multi_read_barcodes = barcode_counts[barcode_counts > 1].index

    # create single and multi read tables
    temp_table_single = full_table[full_table['BC'].isin(single_read_barcodes)].copy()
    temp_table_multi = full_table[full_table['BC'].isin(multi_read_barcodes)].copy()

    # set the max read count, total read count, and mode for single read barcodes
    temp_table_single['mCount'] = 1
    temp_table_single['tCount'] = 1
    temp_table_single['Mode'] = 'Amb'

    return temp_table_single, temp_table_multi


def split_multi_read_barcodes_into_clean_and_chimeric(temp_table_multi: pd.DataFrame)-> tuple:
    """
    Split multi-read barcodes into clean and chimeric based on LUTnr associations.

    :param temp_table_multi: DataFrame containing the multi-read barcodes

    :return: Tuple containing the clean and chimeric multi-read tables
    """
    logger.info("Splitting multi-read barcodes into clean and chimeric")
    temp_table_multi['mismatches'] = temp_table_multi['mismatches'].astype(float)

    # Group by BC and LUTnr to compute statistics
    group_cols = ['BC', 'LUTnr']
    temp_table_multi_grouped = temp_table_multi.groupby(group_cols).agg({
        'bitScore': 'mean',
        'mismatches': 'median',
        'Reads': 'count'
    }).reset_index().rename(columns={'Reads': 'tCount'})

    # Identify clean and chimeric barcodes
    bc_counts = temp_table_multi_grouped['BC'].value_counts()
    clean_barcodes = bc_counts[bc_counts == 1].index
    chimeric_barcodes = bc_counts[bc_counts > 1].index

    temp_table_multi_clean = temp_table_multi_grouped[temp_table_multi_grouped['BC'].isin(clean_barcodes)].copy()
    temp_table_multi_chimeric = temp_table_multi_grouped[temp_table_multi_grouped['BC'].isin(chimeric_barcodes)].copy()

    temp_table_multi_clean['mCount'] = temp_table_multi_clean['tCount']
    temp_table_multi_clean['Mode'] = 'Def'

    return temp_table_multi_clean, temp_table_multi_chimeric

def calculate_consensus_alignment(temp_table_multi_chimeric):
    """
    Calculate consensus alignment of chimeric barcodes.
    """
    logger.info("Calculating consensus alignment for chimeric barcodes")
    # For chimeric barcodes, select the LUTnr with the maximum mCount
    temp_table_multi_chimeric['mCount'] = temp_table_multi_chimeric['tCount']
    temp_table_multi_chimeric['tCount'] = temp_table_multi_chimeric.groupby('BC')['tCount'].transform('sum')

    # Select the LUTnr with the highest mCount for each barcode
    idx = temp_table_multi_chimeric.groupby('BC')['mCount'].idxmax()
    temp_table_multi_consensus = temp_table_multi_chimeric.loc[idx].copy()
    temp_table_multi_consensus['Mode'] = 'Def'

    return temp_table_multi_consensus

def combine_tables(temp_table_multi_final, temp_table_single):
    """
    Combine the multi-read and single-read tables into the final output table.
    """
    logger.info("Combining tables to create final output")
    output_table = pd.concat([temp_table_multi_final, temp_table_single], ignore_index=True)
    return output_table

def main():
    start_time = time.time()

    # load config
    config = get_config(4)

    lut_df = pd.read_csv(config["in_name_LUT"])

    # Load trimmed reads and barcodes
    reads_trim, reads_BC = load_trimmed_reads(config["fragment_file"], config["barcode_file"])

    # Create BLAST database from LUT sequences
    blast_db_prefix = make_customarray_reference_index(lut_df)

    # Save unique fragments as FASTA file
    fragments_unique_fa_name, unique_reads_sequences = save_unique_fragments(reads_trim)

    # Align against the library using BLAST
    blast_out_file_name = align_against_library(fragments_unique_fa_name, blast_db_prefix)

    # Compress BLAST output
    blast_output_gz = f"./0_data/{filename}_blastOutput.csv.gz"
    compress_blast_output(blast_out_file_name, blast_output_gz)

    # Read BLAST output
    blast_df = read_blast_output(blast_out_file_name, unique_reads_sequences, lut_df)

    # Create full table with fragments and barcodes
    full_table = create_full_table(blast_df, lut_df, reads_trim, reads_BC)

    # Perform Starcode barcode reduction
    starcode_exploded = starcode_based_barcode_reduction(barcodes_file)

    # Replace barcodes with Starcode-reduced versions
    full_table = replace_barcodes_with_starcode_versions(full_table, starcode_exploded)

    # Split reads into single-read and multi-read barcodes
    temp_table_single, temp_table_multi = split_reads_into_single_and_multi_read_barcodes(full_table)

    # Split multi-read barcodes into clean and chimeric
    temp_table_multi_clean, temp_table_multi_chimeric = split_multi_read_barcodes_into_clean_and_chimeric(temp_table_multi)

    # Calculate consensus alignment for chimeric barcodes
    temp_table_multi_consensus = calculate_consensus_alignment(temp_table_multi_chimeric)

    # Combine clean and consensus tables
    temp_table_multi_final = pd.concat([temp_table_multi_clean, temp_table_multi_consensus], ignore_index=True)

    # Combine all tables into final output
    output_table = combine_tables(temp_table_multi_final, temp_table_single)

    # Save the output table
    output_file = f"0_data/{filename}_multipleContfragmentsComplete.csv"
    output_table.to_csv(output_file, index=False)
    logger.info(f"Output saved to {output_file}")

    # Print total analysis time
    total_time = time.time() - start_time
    logger.info(f"Total analysis time: {total_time:.2f} seconds")

if __name__ == "__main__":
    main()
