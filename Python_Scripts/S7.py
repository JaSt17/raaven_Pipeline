import os
import pandas as pd
import logging
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
                print(f"Error reading {file}: {e}")

    # Combine the DataFrames if the list is not empty
    if df_list:
        try:
            combined_data = pd.concat(df_list, ignore_index=True)
            return combined_data
        except ValueError as e:
            print(f"Error concatenating DataFrames: {e}")
            return pd.DataFrame()
    else:
        print("No CSV files found or all files failed to read.")
        return pd.DataFrame()
    

def split_reference_name(df: pd.DataFrame) -> pd.DataFrame:
    """
    Splits the 'reference_name' column of a pandas DataFrame into multiple new columns
    and removes unnecessary columns

    Parameters:
        df (pd.DataFrame): The input DataFrame containing a 'reference_name' column.

    Returns:
        pd.DataFrame: The modified DataFrame with new columns added and unnecessary columns removed.
    """
    # Check if 'reference_name' column exists in the DataFrame
    if 'reference_name' not in df.columns:
        raise ValueError("The DataFrame does not have a 'reference_name' column.")

    # Split 'reference_name' into multiple columns
    split_cols = df['reference_name'].str.split(',', expand=True)

    # Expected number of columns after split
    expected_cols = ["Category", "Protein", "Origin", "Extra", "Number", "GeneName"]

    # Handle cases where the number of splits is less than expected
    num_splits = split_cols.shape[1]
    if num_splits < len(expected_cols):
        # Add missing columns with NaN values
        for i in range(num_splits, len(expected_cols)):
            split_cols[i] = None

    # Assign the split columns to new column names
    split_cols.columns = expected_cols

    # Concatenate the new columns to df
    df = pd.concat([df, split_cols], axis=1)

    # Remove unnecessary columns
    columns_to_drop = ['reference_name', 'Protein', 'Origin', 'Extra', 'Number']
    existing_columns_to_drop = [col for col in columns_to_drop if col in df.columns]
    df = df.drop(columns=existing_columns_to_drop)

    return df


def normalize_read_counts(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize the read counts in a DataFrame by adjusting for the total RNA counts in each group.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing read count columns.

    Returns:
        pd.DataFrame: The DataFrame with normalized read counts.
    """
    # create a dictionary to store the sum of RNA counts for each group
    group_RNAcount_dict = dict(df.groupby("Group").sum()["RNAcount"])
    
    # Normalize the group RNA counts by dividing each value by the RNAcounts of the biggest group
    # get the max RNAcount
    max_RNAcount = max(group_RNAcount_dict.values())
    norm_group_RNAcount_dict = {group: RNAcount / max_RNAcount for group, RNAcount in group_RNAcount_dict.items()}
    
    # add the normalized RNA counts to the DataFrame
    df["Normalized_RNAcount"] = df.apply(lambda row: row["RNAcount"] / norm_group_RNAcount_dict[row["Group"]], axis=1)
    
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


def get_subset_include(df: pd.DataFrame, include_substr:list) -> pd.DataFrame:
    """
    Get the subset of the DataFrame after including only groups that contain any of the given string in the list.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing a 'Group' column.
        include_substr (list): A list of substrings to include in the group names.

    Returns:
        pd.DataFrame: The subset of the DataFrame after including only the specified groups.
    """
    # get a list of columns that contain any of the substrings
    groups = df['Group'].unique()
    include_groups = [group for group in groups if any(substring in group for substring in include_substr)]
    
    subset = df[df['Group'].isin(include_groups)]
    
    return subset


def main():
    # Load configuration
    config = get_config("S6")
    
    # Load the combined data from the specified directory
    combined_data = load_combined_data(config["input_dir"])
    if combined_data.empty:
        logger.error("No data loaded. Exiting.")
        return
    
    # Add the library fragment to the combined data
    library_fragments = pd.read_csv(config["library_fragments"])
    library_fragments['Group'] = "DNA_pscAAVlib"
    combined_data = pd.concat([combined_data, library_fragments], ignore_index=True)
    
    # Split the 'reference_name' column into multiple columns
    combined_data = split_reference_name(combined_data)
    
    # Normalize the read counts in the DataFrame
    combined_data = normalize_read_counts(combined_data)
    
    # Get the subset of the DataFrame after excluding certain groups
    DNA_resistant = get_subset_exclude(combined_data, config["exclude_groups1"], "DNAse_resistant")
    mRNA_all = get_subset_exclude(DNA_resistant, config["exclude_groups2"], "mRNA_all")
    
    # combine the two subsets with the combined data
    combined_data = pd.concat([combined_data, DNA_resistant, mRNA_all], ignore_index=True)
    
    
    
    # Save the processed data to a new CSV file
    combined_data.to_csv(config["output_table"], index=False)
    logger.info(f"Saved processed data to {config['output_table']}")
    