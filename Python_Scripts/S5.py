#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script merges the information from the library fragments with the information of the fragments in their reference.

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
    

def main():
    start_time = datetime.now()
    # Load configuration
    config = get_config("S5")
    
    # Create a logger
    create_logger(config["log_dir"], "S5")
    
    # give the in_name_LUT a default value if it is None
    if config["in_name_LUT"] is None:
        config["in_name_LUT"] = "No_file"
    
    # Load the input files and see if there is a reference file
    try:
        library_barcodes = pd.read_csv(config["input_table"], dtype={7: str})
        fragments_pos = pd.read_csv(config["in_name_LUT"], dtype={0: str})
    except FileNotFoundError as e:
        logger.error(f"Could not find input file: {e}")
        # If no reference file was found, save the library fragments with RNAcount
        try:
            library_barcodes = pd.read_csv(config["input_table"], dtype={7: str})
            # add column RNAcount and set it to tcount
            library_barcodes.rename(columns={"tCount": "RNAcount"}, inplace=True)
            library_barcodes.to_csv(config["output_table"], index=False)
            logger.info(f"Could not find fragments position file, saved library fragments with RNAcount to {config['output_table']}")
            sys.exit(0)
        # If the library fragments file was not found, exit the script
        except FileNotFoundError as e:
            logger.error(f"Could not find input file: {e}")
            sys.exit(1)
        sys.exit(1)
    
    logger.info(f"Number of unique fragments in the library: {len(library_barcodes['LUTnr'].unique())}")
    logger.info(f"Number of unique fragments created with the reference file: {len(fragments_pos['LUTnr'].unique())}")
    logger.info(f"Percentage of fragments found in the library: {len(library_barcodes['LUTnr'].unique()) / len(fragments_pos['LUTnr'].unique()) * 100:.2f}%")
    
    # Drop the 'Sequence' column from lut_dna
    fragments_pos.drop(columns=['Sequence'], inplace=True)

    # merge library_barcodes with the LUT 
    library_barcodes = pd.merge(library_barcodes, fragments_pos, on=["LUTnr","Peptide"], how="inner")
    
    # Rename the 'Reads' coulmn to 'Sequence'
    library_barcodes.rename(columns={"Reads": "Sequence"}, inplace=True)
    
    # remove unnecessary columns Reads,identity,alignmentLength,gapOpens,q_start,q_end,s_start,s_end,evalue,
    unessary_columns = ["Reads", "identity", "alignmentLength", "gapOpens", "q_start", "q_end", "s_start", "s_end", "evalue",'bitScore','mismatches',]
    # remove each column if it exists
    for col in unessary_columns:
        if col in library_barcodes.columns:
            library_barcodes.drop(col, axis=1, inplace=True)    
    
    # add RNAcount column
    library_barcodes["RNAcount"] = library_barcodes["tCount"]
    
    # reorder columns
    library_barcodes = library_barcodes[['Origion_seq','Mode','Structure', 'LUTnr',  'BC', 'AAstart', 'AAend', 'Peptide',
                                        'start', 'end', 'width', 'Sequence', 'mCount', 'tCount', 'RNAcount']]
    
    # save the merged library_barcodes
    library_barcodes.to_csv(config["output_table"], index=False)
    logger.info(f"Saved annotated library fragments to {config['output_table']}")
    logger.info(f"Finished processing in {datetime.now() - start_time}")
    
    
if __name__ == "__main__":
    main()