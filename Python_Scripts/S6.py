#!/usr/bin/env python3
"""
Barcoded extraction and reduction from all given Samples.

Author: Jaro Steindorff

This script extracts barcodes from given samples, reduces them using the Starcode algorithm, and identifies
corresponding fragments. It processes a csv file with file_path and basename and saves the results for each base name group.

Workflow:
    - Load necessary data from previous steps (LUT.csv, MatchedFragments.csv, fragment_pos.csv)
    - For each RNA sample:
        - Extract barcodes using bbduk2.sh
        - Reduce barcodes using Starcode
        - Match reduced barcodes with fragments
        - Save found fragments for the sample
    - Save a log table with summary statistics
"""

import sys
import logging
from datetime import datetime
import pandas as pd
# local import
from config import get_config

# Initialize logging with custom format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%H:%M:%S',  # Only show hour, minute, and second
    filemode='w',  # Overwrite log file
    filename='Python_Scripts/Logs/S6.log'  # Log file name
    )
logger = logging.getLogger(__name__)

def main():
    start_time = datetime.now()
    # Load configuration
    config = get_config("S6")
    
    try:
        library_fragments = pd.read_csv(config["input_table"], dtype={7: str})
        fragments_pos = pd.read_csv(config["fragments_pos"])
        LUT = pd.read_csv(config["in_name_LUT"])
    except FileNotFoundError as e:
        logger.error(f"Could not find input file: {e}")
        sys.exit(1)
    
    logger.info(f"Number of unique fragments in the library fragments: {len(library_fragments['LUTnr'].unique())}")
    logger.info(f"Number of unique fragments created with the input file: {len(fragments_pos['LUTnr'].unique())}")
    
    print(len(library_fragments))
    # merge library_fragments with the LUT 
    library_fragments = pd.merge(library_fragments, fragments_pos, on="LUTnr", how="inner")
    print(len(library_fragments))
    
    # remove unnecessary columns Reads,identity,alignmentLength,gapOpens,q_start,q_end,s_start,s_end,evalue,
    unessary_columns = ["Reads", "identity", "alignmentLength", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue"]
    library_fragments = library_fragments.drop(columns=unessary_columns)
    
    # add column RNA_count and set it to tcount
    library_fragments["RNA_count"] = library_fragments["tCount"]
    
    # save the merged library_fragments
    library_fragments.to_csv(config["output_table"], index=False)
    logger.info(f"Saved annotated library fragments to {config['output_table']}")
    logger.info(f"Finished processing in {datetime.now() - start_time}")
    
    
if __name__ == "__main__":
    main()