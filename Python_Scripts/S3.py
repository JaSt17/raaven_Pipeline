#!/usr/bin/env python3
"""
Author: Jaro Steindorff

Inputs for the script are:
    - fragment_seq_file: Path to the file with all fragment sequences
    - constitutive_backbone_sequences: List with the constitutive backbone sequences
    - linker_dict: Dictionary with structure names as keys and linker sequences as values
    - length_dict: Dictionary with structure names as keys and sequence lengths as values
    - trim_dict: Dictionary with structure names as keys and trim ranges as values
    - out_name_LUT: Path to the output file with the LUT data
    - out_name: Path to the output file with the alignment data

Output of the script is:
    - The LUT data with annotated structures
    - The alignment data with sequences
"""

import time
import tempfile
import subprocess
import multiprocessing
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
import pysam
import logging
# local import
from costum_functions import aatodna, load_wSet
from config import get_config

# Initialize logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_LUT(LUT_file: str, backbone_seq: list) -> pd.DataFrame:
    """
    Create a lookup table (LUT) from a file with all fragment sequences.
    
    :param LUT_file: Path to the sequence file
    :return: LUT dataframe
    """
    # Load sequences
    LUT_dna = pd.read_csv(LUT_file, sep="\t", header=0)
    # Remove constitutive backbone sequences
    LUT_dna['Sequence'] = LUT_dna['Sequence'].str.replace(backbone_seq[0], '', regex=False)
    LUT_dna['Sequence'] = LUT_dna['Sequence'].str.replace(backbone_seq[1], '', regex=False)
    # Convert sequences to uppercase
    LUT_dna['Sequence'] = LUT_dna['Sequence'].str.upper()
    # Remove duplicates
    LUT_dna = LUT_dna.drop_duplicates()
    # Create a unique identifier for each sequence
    LUT_dna['LUTnr'] = ['seq_' + str(i+1) for i in range(len(LUT_dna))]
    return LUT_dna


def split_and_annotate_sequences(LUT_dna: pd.DataFrame, linker: dict, length: dict) -> tuple:
    """
    Split sequences based on linker and length, and annotate structures.

    :param LUT_dna: Processed LUT data
    :param linker: Dictionary with structure names as keys and linker sequences as values
    :param length: Dictionary with structure names as keys and sequence lengths as values
    :return: Tuple containing a dictionary of subsets and the annotated LUT_dna DataFrame
    """
    subsets = {}
    # Add 'Structure' column
    LUT_dna['Structure'] = np.nan

    # Initialize a Series to keep track of unprocessed sequences
    unprocessed = pd.Series(True, index=LUT_dna.index)

    # Process sequences based on linker
    for structure, seq_linker in linker.items():
        subset = LUT_dna[unprocessed & LUT_dna['Sequence'].str.startswith(seq_linker)].copy()
        subsets[structure] = subset
        logger.info(f"Number of {structure} sequences (by linker): {len(subset)}")
        # Annotate 'Structure' in LUT_dna
        LUT_dna.loc[subset.index, 'Structure'] = structure
        # Mark sequences as processed
        unprocessed.loc[subset.index] = False

    # Process sequences based on length
    for structure, seq_length in length.items():
        subset = LUT_dna[unprocessed & (LUT_dna['Sequence'].str.len() == seq_length)].copy()
        subsets[structure] = subset
        logger.info(f"Number of {structure} sequences (by length): {len(subset)}")
        # Annotate 'Structure' in LUT_dna
        LUT_dna.loc[subset.index, 'Structure'] = structure
        # Mark sequences as processed
        unprocessed.loc[subset.index] = False

    # Remaining unprocessed sequences
    remaining = LUT_dna[unprocessed]
    subsets['unprocessed'] = remaining
    logger.info(f"Number of unprocessed sequences: {len(remaining)}")

    return subsets, LUT_dna


def trim_sequences(subsets: dict, trim_dict) -> dict:
    """
    Trim sequences based on structure and trim ranges.
    
    :param subsets: Dictionary of subsets
    :param trim_dict: Dictionary with structure names as keys and trim ranges as values
    
    :return: Dictionary of trimmed subsets"""

    for structure, trim_range in trim_dict.items():
        subset = subsets[structure]
        subset['Sequence'] = subset['Sequence'].str.slice(*trim_range)
        subsets[structure] = subset
    
    return subsets


def save_fasta_files(subsets: dict) -> dict:
    """
    Save sequences to FASTA files.

    :param subsets: Dictionary of subsets
    :return: Dictionary of FASTA files
    """

    # Write sequences to FASTA files
    def write_fasta(df, file):
        records = [SeqRecord(Seq(seq), id=lutnr, description='') for seq, lutnr in zip(df['Sequence'], df['LUTnr'])]
        SeqIO.write(records, file.name, "fasta")

    fa_files = {}
    for structure, subset in subsets.items():
        # Create temporary file
        fa_file = tempfile.NamedTemporaryFile(prefix=f"LUT_{structure}_", suffix=".fa", delete=False)
        # Write sequences to FASTA file
        write_fasta(subset, fa_file)
        fa_files[structure] = fa_file

    return fa_files


def build_bowtie_index(LUT_file: str):
    """
    Build Bowtie2 index from original sequences.
    
    :return: Path to Bowtie2 index
    """
    # Read original sequences
    seqs_original = list(SeqIO.parse(LUT_file, "fasta"))
    # Translate sequences to amino acids
    seqs_AA = [seq_record.seq.translate() for seq_record in seqs_original]
    # Reverse translate amino acids to DNA using human codon usage
    ids = [seq_record.description for seq_record in seqs_original]
    seqs_optimized = [aatodna(str(aa_seq)) for aa_seq in seqs_AA]
    seqs_optimized_records = [SeqRecord(Seq(seq), id=seq_id, description='') for seq, seq_id in zip(seqs_optimized, ids)]
    # Write to FASTA file
    bowtie_fasta = tempfile.NamedTemporaryFile(prefix="bowtie_", suffix=".fa", delete=False)
    SeqIO.write(seqs_optimized_records, bowtie_fasta.name, "fasta")
    # Build Bowtie2 index
    bowtie_idx_prefix = tempfile.mktemp(prefix="IDX_bowtie_", dir=tempfile.gettempdir())
    subprocess.run(["bowtie2-build", bowtie_fasta.name, bowtie_idx_prefix], check=True)
    return bowtie_idx_prefix


def align_to_reference(fa_file: str, bowtie_idx_prefix: str, structure_name: str):
    """
    Align sequences to reference using Bowtie2.

    :param fa_file: Path to FASTA file with sequences
    :param bowtie_idx_prefix: Path to Bowtie2 index
    :param structure_name: Name of the structure

    :return: DataFrame with alignment information
    """
    # Create temporary files and names
    name_bowtie = tempfile.mktemp(prefix="bowtie_", dir=tempfile.gettempdir())
    sam_file = name_bowtie + ".sam"
    bam_file = name_bowtie + ".bam"
    sorted_bam_file = name_bowtie + "_sort.bam"
    # Get number of threads
    num_threads = multiprocessing.cpu_count()
    # Run Bowtie2 alignment
    bowtie2_cmd = [
        "bowtie2", "--non-deterministic", "--threads", str(num_threads),
        "--very-sensitive", "-f", "-a",
        "-x", bowtie_idx_prefix, "-U", fa_file, "-S", sam_file
    ]
    # Run Bowtie2
    subprocess.run(bowtie2_cmd, check=True)
    # Convert SAM to BAM and sort
    subprocess.run(["samtools", "view", "-@", str(num_threads), "-Sb", sam_file, "-o", bam_file], check=True)
    subprocess.run(["samtools", "sort", "-@", str(num_threads), bam_file, "-o", sorted_bam_file], check=True)
    # Read sorted BAM file
    bam = pysam.AlignmentFile(sorted_bam_file, "rb")
    # Extract alignment information
    frag_ranges = list(bam.fetch(until_eof=True))
    bam.close()
    # Extract alignment information
    alignment_data = []
    for aln in frag_ranges:
        alignment_data.append({
            'LUTnr': aln.query_name,
            'reference_name': aln.reference_name,
            'flag': aln.flag,
            'reference_start': aln.reference_start,
            'mapping_quality': aln.mapping_quality,
            'cigarstring': aln.cigarstring,
            'is_reverse': aln.is_reverse,
            'structure': structure_name
        })
    return pd.DataFrame(alignment_data)


def main():
    # Start timer
    start_time = time.time()
    # load the configuration
    config = get_config("S3")
    # Load and process LUT data from file
    LUT_dna = create_LUT(config["fragment_seq_file"], config["constitutive_backbone_sequences"])
    # Split sequences based on conditions (linker and length) & annotate structures
    subsets, LUT_dna = split_and_annotate_sequences(LUT_dna, config["linker_dict"], config["length_dict"])
    # Save LUT data
    LUT_dna.to_csv(config["out_name_LUT"], index=False)
    # Trim sequences
    subsets = trim_sequences(subsets, config["trim_dict"])
    # Save sequences to FASTA files
    fa_files_dict = save_fasta_files(subsets)
    # Build Bowtie2 index
    bowtie_idx_prefix = build_bowtie_index()
    # Align sequences and process results
    bowtie_results_dfs = []
    for structure, fa_file in fa_files_dict.items():
        df = align_to_reference(fa_file.name, bowtie_idx_prefix, structure)
        bowtie_results_dfs.append(df)
    # Concatenate all alignment results
    df_all_fragments = pd.concat(bowtie_results_dfs)
    # Merge with LUT_dna to get sequences
    df_all_fragments = df_all_fragments.merge(LUT_dna[['LUTnr', 'Sequence']], on='LUTnr', how='left')
    # Save final data
    df_all_fragments.to_csv(config["out_name"], index=False)
    # Print total execution time
    logger.info("Total execution time:")
    logger.info(time.time() - start_time)

if __name__ == "__main__":
    main()
