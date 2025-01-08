#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script extracts barcodes and fragments from sequencing files, matches them to the reference, pairs them, and saves them in separate files.

Workflow:
    - Use bbduk2.sh to extract barcodes from the given barcode sequencing file
    - Use bbduk2.sh to extract fragments from the given fragment sequencing file
    - Use vsearch to match the fragment reads to the reference sequence
    - Pair the barcodes and fragments using seqkit pair
    - Save the paired barcodes and fragments in separate files

Inputs for the script:
    - in_name_barcode: Path to the barcode sequencing file
    - in_name_fragment: Path to the fragment sequencing file
    - out_dir: Output directory path
    - out_name: Output files name
    
Outputs of the script:
    - barcode_out_name: Paired barcode file
    - fragment_out_name: Paired fragment file
"""

import os
import sys
import shutil
import subprocess
import tempfile
import re
import logging
from datetime import datetime
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
    # ensure that the path exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)
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


def main():
    # Start time
    start_time = datetime.now()
    # Get configuration
    config = get_config("S2")
    # Create logger
    create_logger(config["log_dir"], "S2")

    # Get input parameters
    in_name_barcode, in_name_fragment = config["in_name_barcode"], config["in_name_fragment"]
    out_dir, out_name = config["out_dir"], config["out_name"]
    
    # get the number of available threads
    threads = os.cpu_count()
    
    os.makedirs(out_dir, exist_ok=True)
    out_name_barcode = os.path.join(out_dir, f"barcode_{out_name}.fastq.gz")
    bbduk2_args_BC = config["bbduk2_args_BC"] + [
        f"threads={threads}",
        f"in={in_name_barcode}",
        f"out={out_name_barcode}",
    ]
    _, stderr = run_command(["bbmap/bbduk2.sh"] + bbduk2_args_BC, "bbduk2 barcode extraction")
    
    # Extract summary from the barcodes extraction
    summary = extract_summary(stderr)
    if summary:
        logger.info(f"bbduk2 barcode extraction summary:\n{summary}")

    # BBDuk extraction for fragments
    out_name_fragment = os.path.join(out_dir, f"fragment_{out_name}.fastq.gz")
    bbduk2_args_Frag = config["bbduk2_args_Frag"] + [
        f"threads={threads}",
        f"in={in_name_fragment}",
        f"out={out_name_fragment}",
    ]
    _, stderr = run_command(["bbmap/bbduk2.sh"] + bbduk2_args_Frag, "bbduk2 fragment extraction")
    
    # Extract summary from the fragments extraction
    summary = extract_summary(stderr)
    if summary:
        logger.info(f"bbduk2 fragment extraction summary:\n{summary}")
        
    # Use vsearch to match the fragment reads to the reference sequence
    logger.info("Running vsearch to match fragment reads to the reference sequence.")
    
    # set the reference sequence to the human codon-optimized reference sequence
    db = config["input_file"].replace(".fasta", "_HCO.fasta")
    
    with tempfile.NamedTemporaryFile(delete=True, mode='w+', suffix='.txt') as vsearch_out, \
        tempfile.NamedTemporaryFile(delete=True, mode='w+', suffix='.txt') as keep_fragments:
        
        # Run vsearch and store the output in a temporary file
        vsearch_command = [
            f"zcat {out_name_fragment} | "
            f"vsearch --usearch_global - "
            f"--db {db} "
            "--id 1.0 "
            "--minseqlength 1"
            f"--blast6out {vsearch_out.name} "
            f"--threads {threads}"
        ]
        _, stderr = run_command(vsearch_command, "vsearch", shell=True)

        match = re.search(r"(\d+ of \d+ \(\d+\.\d+%\))", stderr)
        if match:
            info = match.group(0)
            logger.info(f"Ratio of matching fragment reads to the reference sequence: {info}")

        # Use awk to extract the matching fragment reads from the vsearch output
        awk_command = [f"awk '{{print $1}}' {vsearch_out.name} > {keep_fragments.name}"]
        _, stderr = run_command(awk_command, "awk", shell=True)

        # Use seqkit to extract the matching fragment reads
        seqkit_command = [
            f"seqkit grep -f {keep_fragments.name} {out_name_fragment} "
            f"-o filtered_fragments.fastq.gz -j {threads}"
        ]
        _, stderr = run_command(seqkit_command, "seqkit grep", shell=True)

        # Move the filtered fragment reads to the output directory
        shutil.move("filtered_fragments.fastq.gz", out_name_fragment)
    

    # Use seqkit pair
    seqkit_command = [f"seqkit pair -1 {out_name_barcode} -2 {out_name_fragment} -u -j {threads}"]

    logger.info("Running seqkit pair to synchronize paired-end reads.")
    _, stderr = run_command(seqkit_command, "seqkit pair", shell=True)
    
    # Regular expressions
    paired_reads_regex = r"(\d+)\s+paired-end reads"
    unpaired_reads_regex = r"(\d+)\s+unpaired reads"

    # Extract paired reads
    paired_reads_match = re.search(paired_reads_regex, stderr)
    paired_reads = int(paired_reads_match.group(1)) if paired_reads_match else 0

    # Extract unpaired reads with file paths
    unpaired_reads_matches = re.findall(unpaired_reads_regex, stderr)

    # Print results
    logger.info(f"Paired reads: {paired_reads}")
    logger.info(f"Unpaired barcode reads: {unpaired_reads_matches[0]}")
    logger.info(f"Unpaired fragment reads: {unpaired_reads_matches[1]}")
    
    os.makedirs(out_dir, exist_ok=True)
    paired_barcode_out = os.path.join(out_dir, f"barcode_{out_name}.fastq.gz")
    paired_fragment_out = os.path.join(out_dir, f"fragment_{out_name}.fastq.gz")
    
    # output names 
    seqkit_barcode_out = out_name_barcode.replace(".fastq.gz", ".paired.fastq.gz")
    seqkit_fragment_out = out_name_fragment.replace(".fastq.gz", ".paired.fastq.gz")
    seqkit_barcode_unpaired = out_name_barcode.replace(".fastq.gz", ".unpaired.fastq.gz")
    seqkit_fragment_unpaired = out_name_fragment.replace(".fastq.gz", ".unpaired.fastq.gz")
    
    # use shutil to move the files to the output directory
    shutil.move(seqkit_barcode_out, paired_barcode_out)
    shutil.move(seqkit_fragment_out, paired_fragment_out)
    
    #remove the unpaired files
    os.remove(seqkit_barcode_unpaired)
    os.remove(seqkit_fragment_unpaired)

    logger.info(f"Paired reads saved to: {paired_barcode_out}")
    logger.info(f"Paired fragments saved to: {paired_fragment_out}")
    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()