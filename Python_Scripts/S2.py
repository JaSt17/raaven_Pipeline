#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script extracts barcodes and fragments from sequencing files, pairs them, and saves them in separate files.

Workflow:
    - Use bbduk2.sh to extract barcodes from the P5 sequencing file
    - Use bbduk2.sh to extract fragments from the P7 sequencing file
    - Pair the barcodes and fragments using seqkit for exact matching
    - Save the paired barcodes and fragments in separate files

Inputs for the script:
    - in_name_barcode: Path to the barcode sequencing file
    - in_name_fragment: Path to the P7 sequencing file
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

# Function to create a global logger
def create_logger(path: str, name: str) -> None:
    filename = path + name + ".log"
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(message)s",
        datefmt="%H:%M:%S",
        filemode="w",
        filename=filename,
    )
    global logger
    logger = logging.getLogger(name)

def run_command(command: list, description: str) -> tuple:
    logger.info(f"Running {description}")
    with subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True) as process:
        stdout, stderr = process.communicate()
        if process.returncode != 0:
            logger.error(f"Error running {description} with code {process.returncode}")
            logger.error(f"stdout: {stdout}")
            logger.error(f"stderr: {stderr}")
            sys.exit(1)
        return stdout, stderr

def extract_summary(stdout: str) -> str:
    match = re.search(r"(Input:.*?Result:.*?bases\s*\(\d+\.\d+%\))", stdout, re.DOTALL)
    if match:
        return match.group(1)
    return None

def main():
    start_time = datetime.now()
    config = get_config("S2")
    create_logger(config["log_dir"], "S2")

    in_name_barcode, in_name_fragment = config["in_name_barcode"], config["in_name_fragment"]
    out_dir, out_name = config["out_dir"], config["out_name"]
    threads = os.cpu_count()

    # BBDuk extraction for barcodes
    out_name_P5 = tempfile.NamedTemporaryFile(prefix="BC_", suffix=".fastq.gz", delete=False).name
    bbduk2_args_BC = config["bbduk2_args_BC"] + [
        f"threads={threads}",
        f"in={in_name_barcode}",
        f"out={out_name_P5}",
    ]
    _, stderr = run_command(["bbmap/bbduk2.sh"] + bbduk2_args_BC, "bbduk2 barcode extraction")
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
    summary = extract_summary(stderr)
    if summary:
        logger.info(f"bbduk2 fragment extraction summary:\n{summary}")
    in_name_fragment = out_name_P7

    seqkit_command = [
        "seqkit", "pair",
        "--threads", str(threads),
        "-1", in_name_barcode,
        "-2", in_name_fragment,
        "-u"
    ]

    _, stderr = run_command(seqkit_command, "seqkit pair")
        
    #output names from seqkit
    P5_paired_out = in_name_barcode.replace(".fastq.gz", ".paired.fastq.gz")
    P7_paired_out = in_name_fragment.replace(".fastq.gz", ".paired.fastq.gz")
    P5_unpaired_out = in_name_barcode.replace(".fastq.gz", ".unpaired.fastq.gz")
    P7_unpaired_out = in_name_fragment.replace(".fastq.gz", ".unpaired.fastq.gz")
    
    name_dict = {
        P5_paired_out : f"barcode_{out_name}.fastq.gz",
        P7_paired_out : f"fragment_{out_name}.fastq.gz",
        P5_unpaired_out : f"barcode_unpaired_{out_name}.fastq.gz",
        P7_unpaired_out : f"fragment_unpaired_{out_name}.fastq.gz",
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
