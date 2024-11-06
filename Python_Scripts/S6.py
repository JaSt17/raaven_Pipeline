#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script annotates the fragments found in the library with the information of the fragments position file.

Workflow:
    - Loading both input files
    - Merging the library fragments with the fragments position file based on the LUTnr
    - Removing unnecessary columns
    - Adding the RNAcount column
    - Saving the annotated library fragments
    
Input:
    - input_table: str - Path to the library fragments file
    - fragments_pos: str - Path to the fragments position file
    - output_table: str - Path to save the annotated library fragments
    
Output:
    - annotated library fragments: csv - The library fragments with the fragments position information

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
    
    # add column RNAcount and set it to tcount
    library_fragments["RNAcount"] = library_fragments["tCount"]
    
    # save the merged library_fragments
    library_fragments.to_csv(config["output_table"], index=False)
    logger.info(f"Saved annotated library fragments to {config['output_table']}")
    logger.info(f"Finished processing in {datetime.now() - start_time}")
    
    
if __name__ == "__main__":
    main()