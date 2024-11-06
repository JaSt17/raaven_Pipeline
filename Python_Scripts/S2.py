#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script extracts barcodes and fragments from sequencing files, pairs them, and saves them in separate files.

Workflow:
    - Use bbduk2.sh to extract barcodes from the P5 sequencing file
    - Use bbduk2.sh to extract fragments from the P7 sequencing file
    - Pair the barcodes and fragments using pairfq
    - Save the paired barcodes and fragments in separate files

Inputs for the script:
    - in_name_P5: Path to the P5 sequencing file
    - in_name_P7: Path to the P7 sequencing file
    - out_dir: Output directory path
    - out_name: Output files name
    
Outputs of the script:
    - barcode_out_name: Paired barcode file
    - fragment_out_name: Paired fragment file
"""

import os
import sys
import subprocess
import tempfile
import shutil
import re
import logging
from datetime import datetime
from config import get_config


# Initialize logging with custom format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%H:%M:%S',  # Only show hour, minute, and second
    filemode='w',  # Overwrite log file
    filename='Python_Scripts/Logs/S2.log'  # Log file name
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
            sys.exit(1)

        return stdout, stderr
    
def extract_summary(stdout: str) -> str:
    """
    Extracts the summary from the stdout of a command.
    
    :param stdout: stdout of the command
    """
    match = re.search(r"(Input:.*?Result:.*?bases\s*\(\d+\.\d+%\))", stdout, re.DOTALL)
    if match:
        return match.group(1)
    return None

# Main function
def main():
    start_time = datetime.now()
    config = get_config("S2")

    # File paths and threading info
    in_name_P5, in_name_P7 = config['in_name_P5'], config['in_name_P7']
    out_dir, out_name = config['out_dir'], config['out_name']
    threads = os.cpu_count()

    # BBDuk extraction for barcodes
    # create a temporary file for the barcode extraction
    out_name_P5 = tempfile.NamedTemporaryFile(prefix="BC_", suffix=".fastq.gz", delete=False).name
    # extract the configuration for the barcode extraction and add the input and output files and the number of threads
    bbduk2_args_BC = config['bbduk2_args_BC'] + [f"threads={threads}", f"in={in_name_P5}", f"out={out_name_P5}"]
    # run bbduk2 with the configuration
    _, stderr = run_command(["bbmap/bbduk2.sh"] + bbduk2_args_BC, "bbduk2 barcode extraction")
    
    # Extract and log barcode summary
    summary = extract_summary(stderr)
    if summary:
        logger.info(f"bbduk2 barcode extraction summary:\n{summary}")
    # set the input file for the fragment extraction to the output file of the barcode extraction
    in_name_P5 = out_name_P5

    # BBDuk extraction for fragments
    # create a temporary file for the fragment extraction
    out_name_P7 = tempfile.NamedTemporaryFile(prefix="P7_", suffix=".fastq.gz", delete=False).name
    # extract the configuration for the fragment extraction and add the input and output files and the number of threads
    bbduk2_args_Frag = config['bbduk2_args_Frag'] + [f"threads={threads}", f"in={in_name_P7}", f"out={out_name_P7}"]
    # run bbduk2 with the configuration
    _, stderr = run_command(["bbmap/bbduk2.sh"] + bbduk2_args_Frag, "bbduk2 fragment extraction")
    
    # Extract and log fragment summary
    summary = extract_summary(stderr)
    if summary:
        logger.info(f"bbduk2 fragment extraction summary:\n{summary}")
    # set the input file for the pairing to the output file of the fragment extraction
    in_name_P7 = out_name_P7

    # create a dictionary with the temporary output file names for the paired and unpaired barcodes and fragments
    paired_outs = {
        "P5_paired": tempfile.NamedTemporaryFile(prefix="P5_", suffix=".fastq.gz", delete=False).name,
        "P7_paired": tempfile.NamedTemporaryFile(prefix="P7_", suffix=".fastq.gz", delete=False).name,
        "P5_singlet": tempfile.NamedTemporaryFile(prefix="P5_singlet_", suffix=".fastq.gz", delete=False).name,
        "P7_singlet": tempfile.NamedTemporaryFile(prefix="P7_singlet_", suffix=".fastq.gz", delete=False).name
    }
    pairfq_args = [
        "pairfq", "makepairs", "-c", "gzip",
        "-f", in_name_P5, "-r", in_name_P7,
        "-fp", paired_outs["P5_paired"], "-rp", paired_outs["P7_paired"],
        "-fs", paired_outs["P5_singlet"], "-rs", paired_outs["P7_singlet"], "--stats"
    ]
    stdout, _ = run_command(pairfq_args, "pairfq pairing")
    
    # Extract and log pairing summary
    match = re.search(r"(Total paired reads\s*:\s*\d+.*\nTotal unpaired reads\s*:\s*\d+)", stdout)
    if match:
        last_two_lines = re.sub(r"\s*:\s*", ":\t", match.group(1))
        logger.info(f"pairfq pairing summary:\n{last_two_lines}")

    # Save outputs
    os.makedirs(out_dir, exist_ok=True)
    shutil.move(paired_outs["P5_paired"], os.path.join(out_dir, f"barcode_{out_name}.fastq.gz"))
    shutil.move(paired_outs["P7_paired"], os.path.join(out_dir, f"fragment_{out_name}.fastq.gz"))

    # Cleanup
    for file_path in [in_name_P5, in_name_P7] + list(paired_outs.values()):
        try:
            os.remove(file_path)
        except OSError as e:
            # If the file does not exist, continue
            continue
        
    logger.info(f"Finished extraction of barcodes and fragments")
    logger.info(f"Output files saved in {out_dir}")
    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
