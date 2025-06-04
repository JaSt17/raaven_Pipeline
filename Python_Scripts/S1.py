#!/usr/bin/env python3
"""
Author: Jaro Steindorff

This script reads a FASTA file given by the path in the config file and generates all possible 
fragments of a given length from the nucleotide sequences in that give file. It translates the refrence in human optimized
codon usage and saves that  as a new reference file. Than it creates a Table with information of every fragment and writes it to a CSV file.

Workflow:
    - Read the input FASTA file
    - Generate fragments of a given length from the amino acid sequences
    - Convert the amino acid fragments to DNA fragments
    - Add overhangs to the DNA fragments
    - Write the fragments to a file

Inputs for the script are:
    - input_file: The path to the FASTA file
    - structure_dict: A dictionary with the structure as the key and the length and frequency as the values
    - output_csv: The path to the output CSV file
    - output_name: The path to the output file

Output of the script is:
    - A CSV file containing a lookup table with the columns LUTnr, ID, Peptide, AAstart, AAstop, and Sequence
    - A file with the human codon optimized reference sequences
    - A file containing all the generated fragments with the strutcture, start and end position, and the DNA sequence and protein sequence
"""

import os
from datetime import datetime
from Bio import SeqIO
import pandas as pd
from Bio.Seq import translate
import logging
# local import
from costum_functions import aatodna, load_wSet
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


def read_fasta(file_path: str) -> list:
    """
    This function reads a FASTA file with the reference sequences and returns a list information about the sequences.

    Parameters:
        file_path (str): The path to the FASTA file

    Returns:
        list: A list of dictionaries with the keys ID and Peptide
    """
    all_sequences = list(SeqIO.parse(file_path, "fasta"))
    aa_list = []
    for record in all_sequences:
        this_id = record.description
        this_seq = record.seq
        # Translate the DNA sequence to amino acids using the standard genetic code
        this_aa = translate(this_seq, table=1)
        # change all '*' to 'm' (Met) in the amino acid sequence
        this_aa = str(this_aa).replace("*", "M")
        aa_list.append({
                "ID": this_id,
                "Peptide": str(this_aa)
            })
    return aa_list


def make_all_frags(k: int, aa_list: list, structure:str, length: int, frequency: int) -> list:
    """
    This function generates all possible fragments of a given length from a single amino acid sequence.
    
    Parameters:
        k (int): The index of the amino acid sequence in the aa_list
        aa_list (list): A list of dictionaries with the keys ID and Peptide
        structure (str): The structure of the fragments
        length (int): The length of the fragments
        frequency (int): The frequency at which to generate the fragments
        
    Returns:
        list: A list of dictionaries with the keys ID, Peptide, AAstart, AAstop
    """
    frag_list = []
    unique_frag = set()
    this_full_aa = aa_list[k]["Peptide"]
    for l in range(0, (len(this_full_aa) - length) + 1, frequency):
            this_fragment = this_full_aa[l:l + length]
            if this_fragment not in unique_frag:
                if this_fragment[0] != "M":  # Skipping fragments starting with 'M' (Start codon)
                    frag_list.append({
                        "ID": aa_list[k]["ID"],
                        "AAstart": l,
                        "AAend": l + length,
                        "Structure": structure,
                        "Peptide": this_fragment
                    })
                    unique_frag.add(this_fragment)
    return frag_list


def get_unique_fragments(frag_list: list) -> list:
    """
    This function removes duplicate fragments from a list of fragments.

    Parameters:
        frag_list (list): A list of dictionaries with the keys ID, Peptide, AAstart, and AAstop
        
    Returns:
        list: A list of dictionaries with the keys ID, Peptide, AAstart, and AA
    """
    unique_frag_list = []
    seen_fragments = set()
    for frag in frag_list:
        fragment = frag["Peptide"]
        if fragment not in seen_fragments:
            seen_fragments.add(fragment)
            unique_frag_list.append(frag)
    
    return unique_frag_list


def generate_fragments(aa_list: dict, structure: str, length: int, frequency: int = 1) -> list:
    """
    This function generates all possible fragments of a given length from a list of amino acid sequences.

    Parameters:
        aa_list (list): A list of dictionaries with the keys ID and Peptide
        structure (str): The structure of the fragments
        length (int): The length of the fragments
        frequency (int): The frequency at which to generate the fragments
        
    Returns:
        list: A list of dictionaries with the keys ID, Peptide, AAstart, AAstop, Sequence
    """
    frag_list = []

    # Generate fragments sequentially
    for k in range(len(aa_list)):
        frag_list.extend(make_all_frags(k, aa_list, structure, length, frequency))

    # Filter out fragments containing 'X' (unknown amino acid)
    frag_list = [frag for frag in frag_list if 'X' not in frag["Peptide"]]

    # Get unique fragments for each structure
    frag_list = get_unique_fragments(frag_list)

    # Load codon usage table
    wSet = load_wSet()

    # Convert AA sequences to DNA sequences (human codon optimized)
    for frag in frag_list:
        frag["Sequence"] = aatodna(frag["Peptide"], wSet)

    # Sort the fragments by the DNA sequence
    sorted_fragments = sorted(frag_list, key=lambda x: x["Sequence"])

    return sorted_fragments


def create_LUTnr(df: pd.DataFrame) -> pd.DataFrame:
    """
    Create a lookup table (LUT) from a file with all fragment sequences.
    
    Parameters:
        df (pd.DataFrame): A DataFrame with the columns ID, Peptide, AAstart, AAstop, and Sequence
        
    Returns:
        pd.DataFrame: A DataFrame with the columns LUTnr, ID, Peptide, AAstart, AAstop, and Sequence
    """
    LUT = df.copy()
    # Drop duplicates if there are any
    LUT = LUT.drop_duplicates()
    # Create new columns for the LUTnr
    LUT['LUTnr'] = ['seq_' + str(i+1) for i in range(len(LUT))]
    
    # Create new columns for the LUT 
    # DNA sequnce length, start and end
    LUT['start'] = LUT['AAstart']*3
    LUT['end'] = LUT['AAend']*3
    LUT['width'] = LUT['end'] - LUT['start']
    
    # Rename columns
    LUT = LUT.rename(columns={"ID": "Origion_seq"})

    return LUT


# Main script
def main():
    start_time = datetime.now()

    config = get_config("S1")
    
    # Create a logger
    create_logger(config["log_dir"], "S1")

    # Read the input FASTA file
    aa_list = read_fasta(config["input_file"])
    
    # write a human codon optimized reference file
    # create a dataframe with all sequences form the aa_list
    ref_hco = pd.DataFrame(aa_list)
    # Load codon usage table once
    wSet = load_wSet()
    # Prepare arguments for multiprocessing
    ref_hco['Sequence'] = ref_hco['Peptide'].apply(lambda x: aatodna(x, wSet))
    # write the reference file to a file
    output_file = config["input_file"].replace(".fasta", "_HCO.fasta")
    with open(output_file, "w") as f:
        for i, row in ref_hco.iterrows():
            f.write(f">{row['ID']}\n{row['Sequence']}\n")
    
    # Generate fragments for each structure
    sorted_fragments = []
    for structure, info in config["structure_dict"].items():
        # Generate fragments for the structure with the given length and frequency
        temp_sorted_fragments = generate_fragments(aa_list, structure, info["length"], info["freq"])
        logger.info(f"Number of unqiue {structure} fragments: {len(temp_sorted_fragments)}")
        sorted_fragments.extend(temp_sorted_fragments)
    
    # get all fragments
    all_fragments = [frag['Peptide'] for frag in sorted_fragments]
    logger.info(f"Number of created unique Amino Acid fragments: {len(set(set(all_fragments)))}")
    
    # Create a Dataframe with the fragments and sort them by 'Structure', and 'AAstart'
    df = pd.DataFrame(sorted_fragments)
    keys = ["Structure","AAstart"]
    # sort the dataframe by the keys
    df = df.sort_values(by=keys)
    
    # Write the all fragments to a file
    with open(config["output_name"], "w") as f:
        # write all Sequences in the df to the file
        f.write("Sequence\n")
        for seq in df["Sequence"]:
            f.write(seq + "\n")
            
    LUT = create_LUTnr(df)
    
    # if there is a LibID add it to the LUT
    try:
        LUT["LibID"] = config["LibID"]
    except KeyError:
        pass
    
    # Write the LUT to a file
    LUT.to_csv(config["output_csv"], index=False)
    
    logger.info(f"Total execution time: {datetime.now() - start_time}")

if __name__ == "__main__":
    main()