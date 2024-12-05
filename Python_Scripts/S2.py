#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script extracts barcodes and fragments from sequencing files, pairs them, and saves them in separate files.

Workflow:
    - Use bbduk2.sh to extract barcodes from the given barcode sequencing file
    - Use bbduk2.sh to extract fragments from the given fragment sequencing file
    - Pair the barcodes and fragments using fastq_pair
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
import subprocess
import tempfile
import shutil
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
    
    # BBDuk extraction for barcodes
    out_name_P5 = tempfile.NamedTemporaryFile(prefix="BC_", suffix=".fastq.gz", delete=False).name
    bbduk2_args_BC = config["bbduk2_args_BC"] + [
        f"threads={threads}",
        f"in={in_name_barcode}",
        f"out={out_name_P5}",
    ]
    _, stderr = run_command(["bbmap/bbduk2.sh"] + bbduk2_args_BC, "bbduk2 barcode extraction")
    
    # Extract summary from the barcodes extraction
    summary = extract_summary(stderr)
    if summary:
        logger.info(f"bbduk2 barcode extraction summary:\n{summary}")
    in_name_barcode = out_name_P5

    # BBDuk extraction for fragments
    out_name_P7 = tempfile.NamedTemporaryFile(prefix="Frag_", suffix=".fastq.gz", delete=False).name
    bbduk2_args_Frag = config["bbduk2_args_Frag"] + [
        f"threads={threads}",
        f"in={in_name_fragment}",
        f"out={out_name_P7}",
    ]
    _, stderr = run_command(["bbmap/bbduk2.sh"] + bbduk2_args_Frag, "bbduk2 fragment extraction")
    
    # Extract summary from the fragments extraction
    summary = extract_summary(stderr)
    if summary:
        logger.info(f"bbduk2 fragment extraction summary:\n{summary}")
    in_name_fragment = out_name_P7

    # Pair barcodes and fragments using fastq_pair
    fastq_pair_command = [
        "fastq_pair",
        in_name_barcode,
        in_name_fragment
    ]

    _, stderr = run_command(fastq_pair_command, "fastq_pair")

    # Change the names of the output files
    P5_paired_out = in_name_barcode.replace(".fastq.gz", "_paired_1.fastq")
    P7_paired_out = in_name_fragment.replace(".fastq.gz", "_paired_2.fastq")
    P5_unpaired_out = in_name_barcode.replace(".fastq.gz", "_unmatched_1.fastq")
    P7_unpaired_out = in_name_fragment.replace(".fastq.gz", "_unmatched_2.fastq")

    # Update names for consistency
    name_dict = {
        P5_paired_out: f"barcode_{out_name}.fastq",
        P7_paired_out: f"fragment_{out_name}.fastq",
        P5_unpaired_out: f"barcode_unpaired_{out_name}.fastq",
        P7_unpaired_out: f"fragment_unpaired_{out_name}.fastq",
    }


    # Save outputs
    os.makedirs(out_dir, exist_ok=True)
    for name in name_dict:
        if os.path.exists(name):
            shutil.move(name, os.path.join(out_dir, name_dict[name]))

    logger.info(f"Finished extraction of barcodes and fragments")
    logger.info(f"Output files saved in {out_dir}")
    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()
