#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script processes the found fragments and combines them into a single dataset. This final dataet holds infomration for every fragment
including in the diffrent sampeles. The script also normalizes the read counts in the dataset to adjust for the total RNA counts in each group.

Workflow:
    - Combined all output data from the specified directory into one DataFrame
    - Add the library fragments to the combined data
    - add the reference sequence lengths to the combined data
    - Split the 'reference_name' column to 'Category' and 'GeneName'
    - Normalize the read counts in the DataFrame to adjust for the total RNA counts 
    - Combine information of identical fragments in a DataFrame and aggregate the data to get the following infomration
        - bitScore: weighted average of bitScore
        - tCount: sum of tCount
        - mismatches: median of mismatches
        - mCount: sum of mCount
        - BC: unique BCs
        - LUTnr: unique LUTnrs
        - RNAcount: sum of RNAcounts
        - Normalized_RNAcount: sum of Normalized_RNAcounts
        - log2NormCount: log2 of Normalized_RNAcount
    for each unique fragment
    - Save the processed data to a new CSV file
    
Input:
    - Directory containing CSV files with fragment information
    - CSV file with library fragments
    - FASTA file with reference sequences
    - trim_dict: A dictionary with structure names as keys and lists of start and end positions as values

    
Output:
    - CSV file with the processed data
"""

import os
import glob
import re
from datetime import datetime
import pandas as pd
from Bio.Seq import translate
import numpy as np
import logging
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
# local import
from config import get_config

# Initialize logging with custom format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%H:%M:%S',  # Only show hour, minute, and second
    filemode='w',  # Overwrite log file
    filename='Python_Scripts/Logs/S7.log'  # Log file name
    )
logger = logging.getLogger(__name__)


def load_combined_data(dir_path: str) -> pd.DataFrame:
    """
    Reads all CSV files in the specified directory and combines them into a single DataFrame.

    Parameters:
        dir_path (str): The path to the directory containing CSV files.

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
            df = pd.read_csv(file_path)
            df['Group'] = group_name
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
        return combined_data
    except ValueError as e:
        logger.info(f"Error combining data: {e}")
        return pd.DataFrame()
    
    
def get_ref_sequence_length_df(file_path:str) -> pd.DataFrame:
    """
    Reads a FASTA file containing reference sequences and returns a DataFrame with sequence ids and lengths.
    Args:
        file_path (str): path to the file containing the reference sequences

    Returns:
        pd.DataFrame: DataFrame with the reference sequence ids and lengths
    """
    
    # read in file with original sequences
    seqs_original = list(SeqIO.parse(file_path, "fasta"))
    # translate sequences to amino acids
    seq = [str(seq_record.seq).strip() for seq_record in seqs_original]
    # get sequence ids
    ids = [seq_record.description for seq_record in seqs_original]
    # change id to have underscore instead of space
    ids = [seq_id.replace(" ", "_") for seq_id in ids]
    
    # create a DataFrame with sreference sequence ids and lengths
    ref_seq_len_df = pd.DataFrame({"reference_name": ids, "seqlength": [len(s) for s in seq]})
    
    return ref_seq_len_df
    

def split_reference_name(df: pd.DataFrame) -> pd.DataFrame:
    """
    Splits the 'reference_name' column into 'Category' and 'GeneName',
    and removes the 'reference_name' column.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing a 'reference_name' column.

    Returns:
        pd.DataFrame: The modified DataFrame with new columns added and 'reference_name' removed.
    """
    # Ensure 'reference_name' column exists
    if 'reference_name' not in df.columns:
        raise ValueError("The DataFrame does not have a 'reference_name' column.")
    
    # get all unique reference names
    unique_ref_names = df['reference_name'].unique()
    
    # create the Category and GeneName columns for each unique reference name
    Category = []
    GeneName = []
    for ref_name in unique_ref_names:
        # split the reference name into Category and GeneName
        parts = ref_name.split(",")
        Category.append(parts[0])
        GeneName.append(parts[5])
        
    # create a DataFrame with the unique reference names, Category, and GeneName
    ref_name_df = pd.DataFrame({"reference_name": unique_ref_names, "Category": Category, "GeneName": GeneName})
    
    # merge the reference name DataFrame with the original DataFrame based on the reference name
    df = pd.merge(df, ref_name_df, how="left", on="reference_name")

    # Drop 'reference_name'
    df = df.drop(columns=['reference_name'])

    return df


def normalize_read_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize the read counts in a DataFrame by adjusting for the total RNA counts in each group.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing read count columns.

    Returns:
        pd.DataFrame: The DataFrame with normalized read counts.
    """
    # Sum RNA counts per group and get the maximum RNAcount
    group_totals = df.groupby("Group")["RNAcount"].transform("sum")
    max_RNAcount = group_totals.max()
    
    # Calculate normalized counts in a vectorized way
    df["Normalized_RNAcount"] = df["RNAcount"] / (group_totals / max_RNAcount)
    
    return df


def create_subsets(df: pd.DataFrame, subsets: dict) -> pd.DataFrame:
    """
    Create subsets based on the specified conditions.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing fragment information.
        subsets (dict): A dictionary with subset names as keys and conditions as values.

    Returns:
        pd.DataFrame: A modified DataFrame with the subsets added.
    """
    
    for subset_name, conditions in subsets.items():
        logger.info(f"Creating subset '{subset_name}' with conditions: {conditions}")
        if conditions[0] == 'exclude':
            temp_df = df[~df['Group'].isin(conditions[1:])].copy()  # Explicit copy to avoid warning
            temp_df.loc[:, 'Group'] = subset_name
            df = pd.concat([df, temp_df], ignore_index=True)
        elif conditions[0] == 'include':
            temp_df = df[df['Group'].isin(conditions[1:])].copy()  # Explicit copy to avoid warning
            temp_df.loc[:, 'Group'] = subset_name
            df = pd.concat([df, temp_df], ignore_index=True)
        else:
            raise ValueError(f"Invalid condition: {conditions[0]}")
    
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
    # Create a new column for weighted bitScore calculation
    df['bitScore_times_tCount'] = df['bitScore'] * df['tCount']
    
    # Define aggregation functions
    aggregations = {
        'bitScore_times_tCount': 'sum',
        'tCount': 'sum',
        'mismatches': 'median',
        'mCount': 'sum',
        'BC': lambda x: ','.join(x.unique()),
        'LUTnr': lambda x: ','.join(x.unique()),
        'RNAcount': 'sum',
        'Normalized_RNAcount': 'sum',
    }
    
    # Perform groupby aggregation
    combined_data = df.groupby(key_cols).agg(aggregations).reset_index()
    
    # Calculate the weighted average of bitScore
    combined_data['bitScore'] = combined_data['bitScore_times_tCount'] / combined_data['tCount']
    
    # Compute log2NormCount
    combined_data['log2NormCount'] = np.log2(combined_data['Normalized_RNAcount'])
    
    # Drop the intermediate column
    combined_data.drop(columns=['bitScore_times_tCount'], inplace=True)
    
    # Adjust the start, end, width, and seqlength columns
    for col in ['start', 'end', 'width', 'seqlength']:
        combined_data[col] = np.floor(combined_data[col] / 3).astype(int)
    
    # Calculate AA_pos and AA_rel_pos
    combined_data['AA_pos'] = np.floor(combined_data['start'] + combined_data['width'] / 2).astype(int)
    combined_data['AA_rel_pos'] = combined_data['AA_pos'] / combined_data['seqlength']
    
    return combined_data


def cut_overhangs_vectorized(df: pd.DataFrame, trim_dict: dict) -> pd.DataFrame:
    """ 
    This function cuts the overhangs of the sequences based on the structure of the fragment.
    
    Args:
        df (pd.DataFrame): The input DataFrame containing the 'structure' and 'Sequence' columns.
        trim_dict (dict): A dictionary with structure names as keys and lists of start and end positions as values.
        
    Returns:
        pd.DataFrame: The DataFrame with the overhangs of the sequences cut based on the structure of the fragment.
    """
    for key, value in trim_dict.items():
        df.loc[df['structure'] == key, 'Sequence'] = df.loc[df['structure'] == key, 'Sequence'].str.slice(value[0], value[1])

    return df

    
def main():
    start_time = datetime.now()
    # Load configuration
    config = get_config("S7")
    
    # Load the combined data from the specified directory
    combined_data = load_combined_data(config["input_dir"])
    if combined_data.empty:
        logger.info("No data loaded. Exiting.")
        return
    
    # Add the library fragment to the combined data
    library_fragments = pd.read_csv(config["library_fragments"])
    library_fragments['Group'] = "DNA_pscAAVlib"
    combined_data = pd.concat([combined_data, library_fragments], ignore_index=True)
    
    # Only take the subset that have Mode == "Def" which are definitive Barcodes that have more than 1 read
    combined_data = combined_data[combined_data["Mode"] == "Def"]
    
    # Get the reference sequence lengths
    logger.info("Getting reference sequence lengths")
    ref_seq_len_df = get_ref_sequence_length_df(config["original_seq_file"])
    # add the reference sequence lengths to the combined data with the reference_name as the key
    combined_data = pd.merge(combined_data, ref_seq_len_df, how="left", on="reference_name")
    
    logger.info("Splitting reference names and renaming columns")
    # Split the 'reference_name' column into multiple columns
    combined_data = split_reference_name(combined_data)
    
    logger.info("Creating Subsets")
    # Create subsets based on the specified conditions
    combined_data = create_subsets(combined_data, config["subsets"])
    
    logger.info("Normalizing read counts")
    # Normalize the read counts in the DataFrame
    combined_data = normalize_read_counts(combined_data)
    
    logger.info("Combining information of identical fragments")
    # Define the key columns for the groupby operation
    key_cols = ["Group", "Category", "GeneName", "structure", "start", "end", "width", "seqlength", "Sequence"]
    # Combine information of identical fragments in a DataFrame
    combined_data = combine_information_of_identical_fragments(combined_data, key_cols)
    
    logger.info("Cutting overhangs of the sequences based on the structure of the fragment")
    # Cut the overhangs of the sequences based on the structure of the fragment
    combined_data = cut_overhangs_vectorized(combined_data, config["trim_dict"])
    
    logger.info("Translating the sequences to amino acids")
    combined_data["Peptide"] = [translate(seq) for seq in combined_data["Sequence"]]
    
    # Save the processed data to a new CSV file
    combined_data.to_csv(config["output_table"], index=False)
    logger.info(f"Saved processed data to {config['output_table']}")
    logger.info(f"Finished processing in {datetime.now() - start_time}")
    
    
if __name__ == "__main__":
    main()
    