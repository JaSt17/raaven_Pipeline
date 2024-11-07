#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script processes the found fragments and combines them into a single dataset.

Workflow:
    - Combined all output data from the specified directory into one DataFrame
    - Add the library fragments to the combined data
    - add the reference sequence lengths to the combined data
    - Split the 'reference_name' column into multiple columns
    - Normalize the read counts in the DataFrame to adjust for the total RNA counts 
    - Get the subset of the DataFrame after excluding certain groups
    - Combine information of identical fragments in a DataFrame
    - Save the processed data to a new CSV file
    
Input:
    - 
    
Output:
    - 
"""

import os
import pandas as pd
from pandarallel import pandarallel
import numpy as np
import logging
from Bio import SeqIO
# local import
from config import get_config

# Initialize logging with custom format
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(message)s',
    datefmt='%H:%M:%S',  # Only show hour, minute, and second
    #filemode='w',  # Overwrite log file
    #filename='Python_Scripts/Logs/S7.log'  # Log file name
    )
logger = logging.getLogger(__name__)


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


def load_combined_data(dir_path: str) -> pd.DataFrame:
    """
    Reads all CSV files in the specified directory and combines them into a single DataFrame.

    Parameters:
        dir_path (str): The path to the directory containing CSV files.

    Returns:
        pd.DataFrame: A DataFrame containing the combined data from all CSV files.
    """
    # Initialize a list to hold DataFrames
    df_list = []
    
    logger.info(f"Loading data from {dir_path}")

    # Loop through the CSV files in the directory
    for file in os.listdir(dir_path):
        if file.endswith(".csv"):
            file_path = os.path.join(dir_path, file)
            # extract the group name from the file name
            group_name = file.split('.')[1]
            try:
                # Read the CSV file
                df = pd.read_csv(file_path)
                # Add a column with the group name
                df['Group'] = group_name
                # Append the DataFrame to the list
                df_list.append(df)
            except Exception as e:
                logger.info(f"Error reading file {file}: {e}")

    # Combine the DataFrames if the list is not empty
    if df_list:
        try:
            combined_data = pd.concat(df_list, ignore_index=True)
            return combined_data
        except ValueError as e:
            logger.info(f"Error combining data: {e}")
            return pd.DataFrame()
    else:
        logger.info("No CSV files found in the directory.")
        return pd.DataFrame()
    

def split_reference_name(df: pd.DataFrame) -> pd.DataFrame:
    """
    Splits the 'reference_name' column of a pandas DataFrame into multiple new columns
    and removes unnecessary columns.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing a 'reference_name' column.

    Returns:
        pd.DataFrame: The modified DataFrame with new columns added and unnecessary columns removed.
    """
    # Ensure 'reference_name' column exists
    if 'reference_name' not in df.columns:
        raise ValueError("The DataFrame does not have a 'reference_name' column.")

    # Split 'reference_name' into columns
    split_cols = df['reference_name'].str.split(',', expand=True)
    
    # Define expected columns and assign split results, filling in missing columns with None
    expected_cols = ["Category", "Protein", "Origin", "Extra", "Number", "GeneName"]
    for i, col_name in enumerate(expected_cols):
        df[col_name] = split_cols[i] if i < split_cols.shape[1] else None

    # Drop unnecessary columns in a single call
    df.drop(columns=['reference_name', 'Protein', 'Origin', 'Extra', 'Number'], errors='ignore', inplace=True)

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


def get_subset_exclude(df: pd.DataFrame, exclude_groups:list, group:str) -> pd.DataFrame:
    """
    Get the subset of the DataFrame after excluding certain groups.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing a 'Group' column.
        exclude_groups (list): A list of group names to exclude from the subset.

    Returns:
        pd.DataFrame: The subset of the DataFrame after excluding the specified groups.
    """
    # Exclude rows with specified group names
    subset = df[~df['Group'].isin(exclude_groups)]
    subset = subset[subset['Group'] == group]
    return subset


def combine_information_of_identical_fragments(df: pd.DataFrame, key_cols: list) -> pd.DataFrame:
    """
    Combine information of identical fragments in a DataFrame.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing fragment information.
        key_cols (list): A list of column names to identify identical fragments.
    
    Returns:
        pd.DataFrame: The DataFrame with combined information for identical fragments.
    """
    # Define a function to perform the required aggregations for each group
    def aggregate_group(df):
        bitScore = (df['bitScore'] * df['tCount']).sum() / df['tCount'].sum()
        mismatches = df['mismatches'].median()
        mCount = df['mCount'].sum()
        tCount = df['tCount'].sum()
        BC = ','.join(df['BC'].unique())
        Animals = ','.join(df['Sample'].unique())
        LUTnrs = ','.join(df['LUTnr'].unique())
        RNAcount = df['RNAcount'].sum()
        NormCount = df['Normalized_RNAcount'].sum()
        
        return pd.Series({
            'bitScore': bitScore,
            'mismatches': mismatches,
            'mCount': mCount,
            'tCount': tCount,
            'BC': BC,
            'Animals': Animals,
            'LUTnrs': LUTnrs,
            'RNAcount': RNAcount,
            'NormCount': NormCount
        })
    # Initialize pandarallel
    pandarallel.initialize(nb_workers=os.cpu_count())
    
    # Optimize groupby operation with parallel apply if necessary
    combined_data = combined_data.groupby(key_cols, as_index=False, group_keys=False).parallel_apply(aggregate_group)
    
    # Adjust the start and width columns for the combined fragments to match AA readingframe
    combined_data['start'] = np.floor(combined_data['start'] / 3).astype(int)
    combined_data['end'] = np.floor(combined_data['end'] / 3).astype(int)
    combined_data['width'] = np.floor(combined_data['width'] / 3).astype(int)
    combined_data['seqlength'] = np.floor(combined_data['seqlength'] / 3).astype(int)
    
    # calculate the absolute AA position of the fragment
    combined_data['AA_pos'] = combined_data['start'] + 1
    # calculate the relative AA position of the fragment
    combined_data['AA_rel_pos'] = combined_data['AA_pos'] / combined_data['seqlength']
    
    return combined_data


def cut_overhangs_vectorized(df: pd.DataFrame) -> pd.DataFrame:
    """ 
    This function cuts the overhangs of the sequences based on the structure of the fragment.
    
    Args:
        df (pd.DataFrame): The input DataFrame containing the 'structure' and 'Sequence' columns.
        
    Returns:
        pd.DataFrame: The DataFrame with the overhangs of the sequences cut based on the structure of the fragment.
    """
    df.loc[df['structure'] == "14aa", 'Sequence'] = df.loc[df['structure'] == "14aa", 'Sequence'].str.slice(2, 44)
    df.loc[df['structure'] == "14aaG4S", 'Sequence'] = df.loc[df['structure'] == "14aaG4S", 'Sequence'].str.slice(14, 56)
    df.loc[df['structure'] == "14aaA5", 'Sequence'] = df.loc[df['structure'] == "14aaA5", 'Sequence'].str.slice(14, 56)
    df.loc[df['structure'] == "22aa", 'Sequence'] = df.loc[df['structure'] == "22aa", 'Sequence'].str.slice(2, 68)
    return df

    
def main():
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
    
    # Get the reference sequence lengths
    logger.info("Getting reference sequence lengths")
    ref_seq_len_df = get_ref_sequence_length_df(config["original_seq_file"])
    # add the reference sequence lengths to the combined data with the reference_name as the key
    combined_data = pd.merge(combined_data, ref_seq_len_df, how="left", on="reference_name")
    
    logger.info("Splitting reference names and renaming columns")
    # Split the 'reference_name' column into multiple columns
    combined_data = split_reference_name(combined_data)
    
    logger.info("Normalizing read counts")
    # Normalize the read counts in the DataFrame
    combined_data = normalize_read_counts(combined_data)
    
    logger.info("Creating DNA_resistant and mRNA_all subsets")
    # get the subset for DNAse_resistant and mRNA
    DNA_resistant = get_subset_exclude(combined_data, config["exclude_groups1"], "DNAse_resistant")
    # get the subset of mRNA_all
    mRNA_all = get_subset_exclude(DNA_resistant, config["exclude_groups2"], "mRNA_all")
    # combine the two subsets with the combined data
    combined_data = pd.concat([combined_data, DNA_resistant, mRNA_all], ignore_index=True)
    
    logger.info("Combining information of identical fragments")
    # Define the key columns
    key_cols = ["Group", "Category", "GeneName", "structure", "start", "width", "Sequence", "seqlength"]

    # Combine information of identical fragments in a DataFrame
    combined_data = combine_information_of_identical_fragments(combined_data, key_cols)
    
    logger.info("Cutting overhangs of the sequences based on the structure of the fragment")
    # Cut the overhangs of the sequences based on the structure of the fragment
    combined_data = cut_overhangs_vectorized(combined_data)
    
    # Save the processed data to a new CSV file
    combined_data.to_csv(config["output_table"], index=False)
    logger.info(f"Saved processed data to {config['output_table']}")
    
    
if __name__ == "__main__":
    main()
    