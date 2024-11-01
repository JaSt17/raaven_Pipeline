#!/usr/bin/env python3
"""
Mapping of Fragment Sequences to Original Sequences to get the position of the fragment in the original sequence

Author: Jaro Steindorff

Inputs for the script are:
    - original_seq_file: Path to the file with all original sequences
    - wSet: Path to the file with the wSet data
    - fragment_seq_file: Path to the file with all fragment sequences
    - constitutive_backbone_sequences: List with the constitutive backbone sequences
    - linker_dict: Dictionary with structure names as keys and linker sequences as values
    - length_dict: Dictionary with structure names as keys and sequence lengths as values
    - trim_dict: Dictionary with structure names as keys and trim ranges as values
    - out_name_LUT: Path to the output file with the LUT data
    - out_name: Path to the output file with the alignment data

Output of the script is:
    - The LUT data with annotated structures
    - The alignment data with reference_name,strand,width,start,end,LUTnr,structure,Sequence
"""

from datetime import datetime
import sys
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
from costum_functions import aatodna , load_wSet
from config import get_config

# Initialize logging with custom format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%H:%M:%S'
)
logger = logging.getLogger(__name__)


def run_command(command: list, description: str) -> tuple:
    """
    Runs a subprocess command and returns stdout, stderr, and error status.
    
    :param command: Command list to execute
    :param description: Description of the command for logging
    """
    logger.info(f"Running {description}")
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            logger.error(f"Error running {description} with code {process.returncode}")
            logger.error(stderr)
            sys.exit(1)

        return stdout, stderr
    

def create_LUT(LUT_file: str, backbone_seq: list) -> pd.DataFrame:
    """
    Create a lookup table (LUT) from a file with all fragment sequences.
    Cuts off the backbone sequences, sets all sequences to uppercase, and adds a unique identifier.
    
    :param LUT_file: Path to the sequence file
    :return: LUT dataframe
    """
    LUT_dna = pd.read_csv(LUT_file, sep="\t", header=0)
    LUT_dna['Sequence'] = LUT_dna['Sequence'].str.replace(backbone_seq[0], '', regex=False)
    LUT_dna['Sequence'] = LUT_dna['Sequence'].str.replace(backbone_seq[1], '', regex=False)
    LUT_dna['Sequence'] = LUT_dna['Sequence'].str.upper()
    LUT_dna = LUT_dna.drop_duplicates()
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
    LUT_dna['Structure'] = np.nan
    LUT_dna['Structure'] = LUT_dna['Structure'].astype(str)
    unprocessed = pd.Series(True, index=LUT_dna.index)

    for structure, seq_linker in linker.items():
        subset = LUT_dna[unprocessed & LUT_dna['Sequence'].str.startswith(seq_linker)].copy()
        subsets[structure] = subset
        logger.info(f"Number of {structure} sequences (by linker): {len(subset)}")
        LUT_dna.loc[subset.index, 'Structure'] = structure
        unprocessed.loc[subset.index] = False

    for structure, seq_length in length.items():
        subset = LUT_dna[unprocessed & (LUT_dna['Sequence'].str.len() == seq_length)].copy()
        subsets[structure] = subset
        logger.info(f"Number of {structure} sequences (by length): {len(subset)}")
        LUT_dna.loc[subset.index, 'Structure'] = structure
        unprocessed.loc[subset.index] = False

    remaining = LUT_dna[unprocessed]
    if len(remaining) > 0:
        subsets['unprocessed'] = remaining
    logger.info(f"Number of unprocessed sequences: {len(remaining)}")

    return subsets, LUT_dna

def trim_sequences(subsets: dict, trim_dict) -> dict:
    """
    Trim sequences based on structure and trim ranges.
    The trim ranges are given as a dictionary with structure names as keys and trim ranges as values.
    
    :param subsets: Dictionary of subsets
    :param trim_dict: Dictionary with structure names as keys and trim ranges as values
    :return: Dictionary of trimmed subsets"""

    for structure, trim_range in trim_dict.items():
        subset = subsets[structure]
        subset['Sequence'] = subset['Sequence'].str.slice(*trim_range)
        subsets[structure] = subset
    
    return subsets

def create_fasta_files(subsets: dict) -> dict:
    """
    Create FASTA files from subsets. The dictionary keys are the structure names.
    These FASTA files are used to build Bowtie2 indices.

    :param subsets: Dictionary of subsets
    :return: Dictionary of FASTA files
    """

    def write_fasta(df, file):
        records = [SeqRecord(Seq(seq), id=lutnr, description='') for seq, lutnr in zip(df['Sequence'], df['LUTnr'])]
        SeqIO.write(records, file.name, "fasta")

    fa_files = {}
    for structure, subset in subsets.items():
        fa_file = tempfile.NamedTemporaryFile(prefix=f"LUT_{structure}_", suffix=".fa", delete=False)
        write_fasta(subset, fa_file)
        fa_files[structure] = fa_file

    return fa_files

def build_bowtie_index(original_sequences: str, wSet_path: str) -> str:
    """
    Build Bowtie2 index from original sequences.

    :param original_sequences: Path to the file with original sequences
    :param wSet_path: Path to the file with the wSet data
    :return: Path to Bowtie2 index
    """
    # read in file with original sequences
    seqs_original = list(SeqIO.parse(original_sequences, "fasta"))
    # translate sequences to amino acids
    seqs_AA = [seq_record.seq.translate() for seq_record in seqs_original]
    # get sequence ids
    ids = [seq_record.description for seq_record in seqs_original]
    # change id to have underscore instead of space
    ids = [seq_id.replace(" ", "_") for seq_id in ids]
    # load wSet with the hsa codon table
    wSet = load_wSet(wSet_path)
    # translate amino acids back to DNA using the hsa specific condon usage table
    seqs_optimized = [aatodna(str(aa_seq), wSet) for aa_seq in seqs_AA]
    # create SeqRecord objects for the optimized sequences
    seqs_optimized_records = [SeqRecord(Seq(seq), id=seq_id, description='') for seq, seq_id in zip(seqs_optimized, ids)]
    # create a temporary FASTA file with the optimized sequences
    bowtie_fasta = tempfile.NamedTemporaryFile(prefix="bowtie_", suffix=".fa", delete=False)
    SeqIO.write(seqs_optimized_records, bowtie_fasta.name, "fasta")
    # create a temporary Bowtie2 index
    bowtie_idx_prefix = tempfile.mktemp(prefix="IDX_bowtie_", dir=tempfile.gettempdir())
    run_command(["bowtie2-build", bowtie_fasta.name, bowtie_idx_prefix], "Bowtie2 index build")

    return bowtie_idx_prefix


def align_to_reference(fa_file: str, bowtie_idx_prefix: str, structure_name: str):
    """
    Align sequences to reference using Bowtie2. The reference is the original sequence file. 
    The alignment information is stored in a DataFrame. SO we get the position of the sequence in the reference.

    :param fa_file: Path to FASTA file with sequences
    :param bowtie_idx_prefix: Path to Bowtie2 index
    :param structure_name: Name of the structure
    :return: DataFrame with alignment information
    """
    name_bowtie = tempfile.mktemp(prefix="bowtie_", dir=tempfile.gettempdir())
    sam_file = name_bowtie + ".sam"
    bam_file = name_bowtie + ".bam"
    sorted_bam_file = name_bowtie + "_sort.bam"
    num_threads = multiprocessing.cpu_count()
    bowtie2_cmd = [
        "bowtie2", "--non-deterministic", "--threads", str(num_threads),
        "--very-sensitive", "-f", "-a",
        "-x", bowtie_idx_prefix, "-U", fa_file, "-S", sam_file
    ]
    run_command(bowtie2_cmd, "Bowtie2 alignment for " + structure_name)
    # check if SAM file is empty
    with open(sam_file, "r") as f:
        if not f.readline():
            logger.warning("SAM file is empty")
            sys.exit(0)
    run_command(["samtools", "view", "-@", str(num_threads), "-Sb", sam_file, "-o", bam_file], "SAM to BAM conversion for " + structure_name)
    run_command(["samtools", "sort", "-@", str(num_threads), bam_file, "-o", sorted_bam_file], "BAM sorting for " + structure_name)
    bam = pysam.AlignmentFile(sorted_bam_file, "rb")
    frag_ranges = list(bam.fetch(until_eof=True))
    bam.close()
    alignment_data = []
    for aln in frag_ranges:
        alignment_data.append({
            'reference_name': aln.reference_name,
            'strand': '-' if aln.is_reverse else '+',
            'width': aln.reference_end - aln.reference_start,
            'start': aln.reference_start,
            'end': aln.reference_end,
            'LUTnr': aln.query_name,
            'structure': structure_name
        })
    return pd.DataFrame(alignment_data)

def main():
    start_time = datetime.now()
    config = get_config("S3")
    LUT_dna = create_LUT(config["fragment_seq_file"], config["constitutive_backbone_sequences"])
    subsets, LUT_dna = split_and_annotate_sequences(LUT_dna, config["linker_dict"], config["length_dict"])
    LUT_dna.to_csv(config["out_name_LUT"], index=False)
    subsets = trim_sequences(subsets, config["trim_dict"])
    fa_files_dict = create_fasta_files(subsets)
    bowtie_idx_prefix = build_bowtie_index(config["original_seq_file"], config["wSet"])
    bowtie_results_dfs = []
    for structure, fa_file in fa_files_dict.items():
        df = align_to_reference(fa_file.name, bowtie_idx_prefix, structure)
        bowtie_results_dfs.append(df)
    df_all_fragments = pd.concat(bowtie_results_dfs)
    df_all_fragments = df_all_fragments.merge(LUT_dna[['LUTnr', 'Sequence']], on='LUTnr', how='left')
    df_all_fragments.to_csv(config["out_name"], index=False)
    
    logger.info(f"LUT data written to {config['out_name_LUT']}")
    logger.info(f"Fragment position data written to {config['out_name']}")
    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
