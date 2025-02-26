#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script processes the found fragments and combines them into a single dataset. This final dataset holds infomration for every found fragment
including in the diffrent samples. The script also normalizes the read counts in the dataset to adjust for the total RNA counts in each group.

Workflow:
    - Combined all output data from the specified directory into one DataFrame
    - Add the library fragments to the combined data
    - add the reference sequence lengths to the combined data
    - Normalize the read counts in the DataFrame to adjust for difference in read depth 
    - Combine information of identical fragments in each group in a DataFrame and aggregate the data to get the following infomration
        - tCount: sum of tCount
        - mCount: sum of mCount
        - BC: unique BCs
        - LUTnr: unique LUTnrs
        - RNAcount: sum of RNAcounts
        - Normalized_RNAcount: sum of Normalized_RNAcounts
        - BC_adjusted_count: log2 of BC_count * Normalized_RNAcount
    - Finally, cut the overhangs of the sequences based on the structure of the fragment (if necessary) and save the processed data to a new CSV file.
    
Input:
    - Directory containing CSV files with fragment information
    - CSV file with library fragments
    - FASTA file with reference sequences
    - trim_dict: A dictionary with structure names as keys and lists of start and end positions as values
    - backbone_seq: A list containing the backbone sequence
    - sample_inputs: The path to the file containing the sample inputs and the corresponding group names
    - output_table: Path to save the processed data
    
Output:
    - CSV file with the processed data
"""

import os
import glob
import re
from datetime import datetime
import pandas as pd
import logging
from concurrent.futures import ThreadPoolExecutor
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


def load_combined_data(dir_path: str, sample_inputs_path: str) -> pd.DataFrame:
    """
    Reads all CSV files in the specified directory and combines them into a single DataFrame.
    Then rename the samples from different organisms to the same name.

    Parameters:
        dir_path (str): The path to the directory containing CSV files.
        sample_inputs_path (str): The path to the file containing the sample inputs and the corresponding group names.

    Returns:
        pd.DataFrame: A DataFrame containing the combined data from all CSV files.
    """
    logger.info(f"Loading data from {dir_path}")
    
    # Use glob to get all CSV files in the directory
    csv_files = glob.glob(os.path.join(dir_path, '*.csv'))
    
    if not csv_files:
        logger.info("No CSV files found in the directory.")
        return pd.DataFrame()

    def read_file(file_path):
        file_name = os.path.basename(file_path)
        # Robust group name extraction using regex
        match = re.match(r'.*?\.(.*?)\..*', file_name)
        group_name = match.group(1) if match else 'Unknown'
        try:
            df = pd.read_csv(file_path, dtype={7: 'str'})  # Example: treat column 7 as strings
            df['Group'] = group_name
            if df.empty:
                logger.info(f"Empty DataFrame in file {file_name}")
                return None
            return df
        except Exception as e:
            logger.info(f"Error reading file {file_name}: {e}")
            return None

    # Use ThreadPoolExecutor
    with ThreadPoolExecutor() as executor:
        # Map the read_file function to the csv_files list
        data_frames = list(executor.map(read_file, csv_files))

    # Filter out None results from failed reads
    data_frames = [df for df in data_frames if df is not None]

    try:
        combined_data = pd.concat(data_frames, ignore_index=True)
    except ValueError as e:
        logger.info(f"Error combining data: {e}")
        return pd.DataFrame()
    
    # load in the sample inputs and the corresponding group names and transform it to a dictionary
    sample_inputs = pd.read_csv(sample_inputs_path)
    # remove the file extension from the sample names
    sample_inputs['Sample'] = sample_inputs['Sample'].str.split(".").str[0]
    sample_inputs_dict = dict(zip(sample_inputs["Sample"], sample_inputs["Group"]))
    
    # rename the samples from different organisms to the same name
    combined_data['Group'] = combined_data['Group'].replace(sample_inputs_dict)
    
    return combined_data
    

def create_subsets(df: pd.DataFrame, subsets: dict) -> pd.DataFrame:
    """
    Create subsets based on the specified conditions.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing fragment information.
        subsets (dict): A dictionary with subset names as keys and conditions as values.
            possible conditions:
                - ('exclude', group1, group2, ...): exclude the specified groups
                - ('include', group1, group2, ...): include only the specified groups
                - ('contains_exclude', string): exclude groups containing the specified string
                - ('contains_include', string): include only groups containing the specified

    Returns:
        pd.DataFrame: A modified DataFrame with the subsets added.
    """
    df_list = [df]
    for subset_name, conditions in subsets.items():
        if conditions[0] == 'exclude':
            temp_df = df[~df['Group'].isin(conditions[1:])].copy()
            temp_df.loc[:, 'Group'] = subset_name
            df_list.append(temp_df)
        elif conditions[0] == 'include':
            temp_df = df[df['Group'].isin(conditions[1:])].copy()  
            temp_df.loc[:, 'Group'] = subset_name
            df_list.append(temp_df)
        elif conditions[0] == 'contains_exclude':
            temp_df = df[~df['Group'].str.contains(conditions[1])].copy()
            temp_df.loc[:, 'Group'] = subset_name
            df_list.append(temp_df)
        elif conditions[0] == 'contains_include':
            temp_df = df[df['Group'].str.contains(conditions[1])].copy()
            temp_df.loc[:, 'Group'] = subset_name
            df_list.append(temp_df)
        else:
            raise ValueError(f"Invalid condition: {conditions[0]}")
    # concatenate the DataFrames
    df = pd.concat(df_list, ignore_index=True)
    
    return df
    

def normalize_read_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize the read counts in the data frame to adjust for difference in read depth.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing read count columns.

    Returns:
        pd.DataFrame: The DataFrame with normalized read counts.
    """
    # Ensure data types are float once
    df = df.astype({'RNAcount': 'float32'})  # Convert 'RNAcount' column to float32 for memory efficiency
    
    # Compute the group-wise sum once and normalize
    group_sum = df.groupby('Group')['RNAcount'].transform('sum')
    df['RNAcount_ratio'] = df['RNAcount'] / group_sum

    return df


def combine_information_of_identical_fragments(df: pd.DataFrame, key_cols: list) -> pd.DataFrame:
    """
    Combine information of identical fragments in a DataFrame.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing fragment information.
        key_cols (list): A list of column names to identify identical fragments.
    
    Returns:
        pd.DataFrame: The DataFrame with combined information for identical fragments.
    """
    
    # Define aggregation functions
    aggregations = {
        'mCount': 'sum',
        'Reads': 'first',
        'BC': lambda x: ','.join(pd.unique(x)),
        'RNAcount': 'sum',
        'RNAcount_ratio': 'sum',
    }
    
    # Perform groupby aggregation
    combined_data = df.groupby(key_cols).agg(aggregations).reset_index()
    
    # Compute the BC_count and BC_ratio
    combined_data['BC_count'] = combined_data['BC'].str.count(',') + 1
    group_bc_count = combined_data.groupby('Group', observed=True)['BC_count'].transform('sum')
    combined_data['BC_ratio'] = combined_data['BC_count'] / group_bc_count

    # Barcode adjusted count ratio
    combined_data['BC_adjusted_count_ratio'] = combined_data['RNAcount_ratio'] + combined_data['BC_ratio'] / 2
    
    # Rename Reads to Sequence
    combined_data.rename(columns={'Reads': 'Sequence'}, inplace=True)

    return combined_data

    
def main():
    start_time = datetime.now()
    # Load configuration
    config = get_config("S6")
    
    # Create a logger
    create_logger(config["log_dir"], "S6")
    
    # Load the combined data from the specified directory
    combined_data = load_combined_data(config["input_dir"], config["sample_inputs"])
    if combined_data.empty:
        logger.info("No data loaded. Exiting.")
        return
    
    # Add the library fragment to the combined data
    library_fragments = pd.read_csv(config["library_fragments"])
    library_fragments['Group'] = config["library_name"]
    # set the columns of the two DataFrames to be the same
    combined_data = combined_data[library_fragments.columns.tolist()]
    # Add the library fragments to the combined data
    combined_data = pd.concat([library_fragments, combined_data], ignore_index=True)
    
    logger.info("Creating Subsets")
    # Create subsets based on the specified conditions
    combined_data = create_subsets(combined_data, config["subsets"])
    
    logger.info("Normalizing read counts")
    # Normalize the read counts in the DataFrame
    combined_data = normalize_read_counts(combined_data)
    
    logger.info("Combining fragment information")
    # Define the key columns for the groupby operation
    key_cols = ["Group", "LUTnr", "Peptide"]
    # Combine information of identical fragments in a DataFrame
    combined_data = combine_information_of_identical_fragments(combined_data, key_cols)
    
    # Save the processed data to a new CSV file
    combined_data.to_csv(config["output_table"], index=False)
    logger.info(f"Saved processed data to {config['output_table']}")
    logger.info(f"Finished processing in {datetime.now() - start_time}")
    
    
if __name__ == "__main__":
    main()
    