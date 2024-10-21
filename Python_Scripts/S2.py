#!/usr/bin/env python3
"""
Extraction of Barcodes and Gene Fragments from the given P5 and P7 raw reads sequencing files

Author: Jaro Steindorff

This script extracts barcodes and fragments from sequencing files, pairs them and saves them in separate files.

Inputs for the script are:
    - in_name_P5: The path to the P5 sequencing file
    - in_name_P7: The path to the P7 sequencing file
    - out_dir: The path to the output directory
    - out_name: The name of the output files
    - bbduk2_args_BC: The arguments for barcode extraction with bbduk2.sh
    - bbduk2_args_Frag: The arguments for fragment extraction with bbduk2.sh

Output of the script is:
    - Two files containing the extracted barcodes and fragments respectively
"""

import os
import sys
import subprocess
import multiprocessing
import tempfile
import shutil
import shlex
import re
import logging
from datetime import datetime
# local import
from config import get_config

# Initialize logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def count_fastq_reads(file_path: str) -> int:
    """
    Count the number of reads in a fastq file.
    
    :param file_path: Path to the fastq file
    :return: Number of reads in the file
    """
    # Command to count the number of sequences
    command = f"pigz -dc {shlex.quote(file_path)} | wc -l| awk '{{print $1/4}}'"

    
    # Run the command
    try:
        result = subprocess.check_output(command, shell=True, text=True)
        # Strip the result and convert it to a float (in case of decimal places)
        num_sequences = int(result.strip())
        return num_sequences
    except subprocess.CalledProcessError as e:
        logger.info(f"Error executing the command: {e}")
        return None

def run_bbduk2(args_list: list) -> tuple:
    """
    Run the bbduk2.sh script with the provided arguments.
    
    :param args_list: List of arguments to pass to bbduk2.sh
    """
    logger.info("Running bbduk2.sh")
    bbduk2_path = os.path.expanduser("~/bbmap/bbduk2.sh")
    # Stream output directly to the console instead of collecting it in memory
    with subprocess.Popen(
        [bbduk2_path] + args_list,
        stdout=sys.stdout,
        stderr=sys.stderr,
        text=True
    ) as process:
        
        process.communicate()

        if process.returncode != 0:
            logger.info(f"Error running bbduk2.sh with return code {process.returncode}", file=sys.stderr)
            sys.exit(1)

        logger.info("Finished bbduk2.sh")

        return process.returncode, process.stderr


def run_pairfq(args_list):
    """
    Run the pairfq tool with the provided arguments.
    
    :param args_list: List of arguments to pass to pairfq"""
    with subprocess.Popen(
        ["pairfq"] + args_list,
        stdout=sys.stdout,
        stderr=sys.stderr,
        text=True
    ) as process:
        process.communicate()

        if process.returncode != 0:
            logger.info(f"Error running pairfq with return code {process.returncode}", file=sys.stderr)
            sys.exit(1)

        logger.info("Finished pairfq")

# Main function
def main():
    # Start timer
    start_time = datetime.now()

    # Read configuration dictionary from the config.py file
    config = get_config("S2")
    in_name_P5 = config['in_name_P5']
    in_name_P7 = config['in_name_P7']
    out_dir = config['out_dir']
    out_name = config['out_name']

    # Count reads
    output_reads = count_fastq_reads(in_name_P5)
    logger.info(f"Utilized Reads: {output_reads} for barcode extraction")

    # Extraction of barcodes
    out_name_P5 = tempfile.NamedTemporaryFile(prefix="BC_", suffix=".fastq.gz", delete=False).name
    threads = multiprocessing.cpu_count()
    # get the bbduk2 arguments for barcode extraction from config
    bbduk2_args_BC = config['bbduk2_args_BC']
    # add the input and output file names to the bbduk2 arguments
    bbduk2_args_BC.extend([
        f"threads={threads}",
        f"in={in_name_P5}",
        f"out={out_name_P5}",
    ])
    logger.info("Running bbduk2 for barcode extraction")
    _, stderr_output = run_bbduk2(bbduk2_args_BC)
    # Parse stderr_output to find the number of reads used for barcode extraction
    match = re.search(r'Result:.*\t(\d+)\sreads', stderr_output)
    if match:
        barcode_read_count = int(match.group(1))
        logger.info(f"Total number of found barcode reads: {barcode_read_count}")

    # Update input file for barcodes
    in_name_P5 = out_name_P5

    # Extraction of fragments
    out_name_P7 = tempfile.NamedTemporaryFile(prefix="P7_", suffix=".fastq.gz", delete=False).name
    # get the bbduk2 arguments for fragment extraction from config
    bbduk2_args_Frag = config['bbduk2_args_Frag']
    bbduk2_args_Frag.extend([
        f"threads={threads}",
        f"in={in_name_P7}",
        f"out={out_name_P7}",
    ])
    logger.info("Running bbduk2 for fragment extraction")
    _, stderr_output = run_bbduk2(bbduk2_args_Frag)
    # Parse stderr_output to find the number of reads used for fragment extraction
    match = re.search(r'Result:.*\t(\d+)\sreads', stderr_output)
    if match:
        fragment_read_count = int(match.group(1))
        logger.info(f"Total number of found fragment reads: {fragment_read_count}")

    # Update input file for fragments
    in_name_P7 = out_name_P7

    # Pairing barcodes and fragments
    out_name_P5_paired = tempfile.NamedTemporaryFile(prefix="P5_", suffix=".fastq.gz", delete=False).name
    out_name_P7_paired = tempfile.NamedTemporaryFile(prefix="P7_", suffix=".fastq.gz", delete=False).name
    out_name_P5_singlet = tempfile.NamedTemporaryFile(prefix="P5_singlet_", suffix=".fastq.gz", delete=False).name
    out_name_P7_singlet = tempfile.NamedTemporaryFile(prefix="P7_singlet_", suffix=".fastq.gz", delete=False).name

    # define the arguments for pairfq
    pairfq_args = [
        "makepairs",
        "-c", "gzip",
        "-f", in_name_P5,
        "-r", in_name_P7,
        "-fp", out_name_P5_paired,
        "-rp", out_name_P7_paired,
        "-fs", out_name_P5_singlet,
        "-rs", out_name_P7_singlet,
        "--stats",
    ]
    logger.info("Running pairfq to pair barcodes and fragments")
    run_pairfq(pairfq_args)

    # Save outputs
    os.makedirs(out_dir, exist_ok=True)
    shutil.move(out_name_P5_paired, f"./{out_dir}/barcode_{out_name}.fastq.gz")
    shutil.move(out_name_P7_paired, f"./{out_dir}/fragment_{out_name}.fastq.gz")

    # Clean up temporary files
    temp_files = [
        in_name_P5,
        in_name_P7,
        out_name_P5_singlet,
        out_name_P7_singlet,
        # Add other temporary files if necessary
    ]
    for temp_file in temp_files:
        try:
            os.remove(temp_file)
        except OSError:
            pass  # File may have been already removed

    # Print total execution time
    logger.info("Total execution time:")
    logger.info(datetime.now() - start_time)

if __name__ == "__main__":
    main()
