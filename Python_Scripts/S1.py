#!/usr/bin/env python3
"""
Costum fragment generation script from input FASTA file

Author: Jaro Steindorff

This script reads a FASTA file given by the path in the config file and generates all possible 
fragments of a given length from the nucleotide sequences in that give file.

Inputs for the script are:
    - input_file: The path to the FASTA file
    - wSet: The path to the codon usage table
    - 14_aa_overhangs: The overhangs to be added to the 14aa fragments
    - 14_aa_G4S_overhangs: The overhangs to be added to the 14aa G4S fragments
    - 14_aa_A5_overhangs: The overhangs to be added to the 14aa A5 fragments
    - 22_aa_overhangs: The overhangs to be added to the 22aa fragments
    - output_name: The path to the output file

Output of the script is:
    - A file containing all the generated fragments
"""

from datetime import datetime
from Bio import SeqIO
from Bio.Seq import translate
from multiprocessing import Pool, cpu_count
import logging
# local import
from costum_functions import aatodna, load_wSet
from config import get_config
import copy


def read_fasta(file_path: str) -> list:
    """
    This function reads a FASTA file and returns a list of dictionaries with the following keys
    Class: The class of the sequence
    Family: The family of the sequence
    Strain: The strain of the sequence
    Note: Any additional notes about the sequence
    Number: The number of the sequence
    Name: The name of the sequence
    AAfragment: The amino acid sequence of the sequence

    :param file_path: The path to the FASTA file

    :return: A list of dictionaries with the keys Class, Family, Strain, Note, Number, Name, and AAfragment
    """
    all_sequences = list(SeqIO.parse(file_path, "fasta"))
    aa_list = []
    for record in all_sequences:
        this_id = record.description
        this_seq = record.seq
        # Translate the DNA sequence to amino acids using the standard genetic code
        this_aa = translate(this_seq, table=1)
        this_id_split = this_id.split(sep=",")
        aa_list.append({
            "Class": this_id_split[0],
            "Family": this_id_split[1],
            "Strain": this_id_split[2],
            "Note": this_id_split[3],
            "Number": this_id_split[4],
            "Name": this_id_split[5],
            "AAfragment": str(this_aa)
        })
    return aa_list


def make_all_frags(k: int, aa_list: list, length: int, frequency: int) -> list:
    """
    This function generates all possible fragments of a given length from a single amino acid sequence.
    
    :param k: The index of the amino acid sequence in the list
    :param aa_list: A list of dictionaries with the keys Class, Family, Strain, Note, Number, Name, and AAfragment
    :param length: The length of the fragments
    :param frequency: The frequency at which to generate the fragments
    
    :return: A list of dictionaries with the keys Class, Family, Strain, Note, Number, Name, AAstart, AAstop, and AAfragment
    """
    frag_list = []
    this_full_aa = aa_list[k]["AAfragment"]
    for l in range(0, (len(this_full_aa) - length) + 1, frequency):
            this_fragment = this_full_aa[l:l + length]
            if this_fragment[0] != "M":  # Skipping fragments starting with 'M' (Start codon)
                frag_list.append({
                    "Class": aa_list[k]["Class"],
                    "Family": aa_list[k]["Family"],
                    "Strain": aa_list[k]["Strain"],
                    "Note": aa_list[k]["Note"],
                    "Number": aa_list[k]["Number"],
                    "Name": aa_list[k]["Name"],
                    "AAstart": l + 1,
                    "AAstop": l + length,
                    "AAfragment": this_fragment
                })
    return frag_list


def get_unique_fragments(frag_list: list) -> list:
    """
    This function removes duplicate fragments from a list of fragments.

    Args:
        frag_list (list): A list of dictionaries with the keys Class, Family, Strain, Note, Number, Name, AAstart, AAstop, and AAfragment

    Returns:
        list: A list of dictionaries with the keys Class, Family, Strain, Note, Number, Name, AAstart, AAstop, and AAfragment
    """
    unique_frag_list = []
    seen_fragments = set()

    for frag in frag_list:
        fragment = frag["AAfragment"]
        if fragment not in seen_fragments:
            seen_fragments.add(fragment)
            unique_frag_list.append(frag)

    return unique_frag_list


def generate_fragments(wSet: str, aa_list: dict, length: int, frequency:int = 1) -> list:
    """
    This function generates all possible fragments of a given length from a list of amino acid sequences.

    :param aa_list: A list of dictionaries with the keys Class, Family, Strain, Note, Number, Name, and AAfragment
    :param min_length: The minimum length of the fragments
    :param max_length: The maximum length of the fragments
    :param frequency: The frequency at which to generate the fragments

    :return: A list of dictionaries with the keys Class, Family, Strain, Note, Number, Name, AAstart, AAstop, and AAfragment
    """

    frag_list = []

    # Use multiprocessing with Pool
    with Pool(cpu_count()) as p:
        # Map the function across all entries in aa_list
        results = p.starmap(make_all_frags, [(k, aa_list, length, frequency) for k in range(len(aa_list))])
    
    # Collect the results
    for result in results:
        frag_list.extend(result)

    # Filter out fragments containing 'X'
    frag_list = [frag for frag in frag_list if not any(c in frag["AAfragment"] for c in "X")]

    # get unique fragments
    frag_list = get_unique_fragments(frag_list)

    # Load codon usage table once
    wSet = load_wSet(wSet)

    # Prepare arguments for multiprocessing
    args_list = [(frag, wSet) for frag in frag_list]

    # Convert AA sequences to DNA sequences using multiprocessing and the 'aatodna' function for human codon optimization
    with Pool(cpu_count()) as p:
        frag_list = p.map(aatodna_parallel, args_list)
    
    # Sort the fragments by 'AAfragment'
    sorted_fragments = sorted(frag_list, key=lambda x: x["DNAfragment"])
    
    return sorted_fragments


def aatodna_parallel(args):
    """
    Wrapper function to convert an amino acid fragment to a DNA fragment using multiprocessing.
    """
    frag, wSet = args
    # Convert AA to DNA
    frag["DNAfragment"] = aatodna(frag["AAfragment"], wSet)
    return frag


# Function to add overhangs to the fragments
def add_overhangs(fragments, five_prime, three_prime):
    """
    This function adds overhangs to the DNA fragments.
    
    :param fragments: A list of dictionaries with the keys Class, Family, Strain, Note, Number, Name, AAstart, AAstop, AAfragment, and DNAfragment
    :param five_prime: The sequence of the 5' overhang
    :param three_prime: The sequence of the 3' overhang
    
    :return: A list of dictionaries with the keys Class, Family, Strain, Note, Number, Name, AAstart, AAstop, AAfragment, and DNAfragment
    """
    for frag in fragments:
        frag["DNAfragment"] = f"{five_prime.lower()}{frag['DNAfragment']}{three_prime.lower()}"
    return fragments


# Main script
def main():
    start_time = datetime.now()

    # Initialize logging with custom format
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(message)s',
        datefmt='%H:%M:%S',  # Only show hour, minute, and second
        filemode='w',  # Overwrite log file
        filename='Python_Scripts/Logs/S1.log'  # Log file name
    )
    logger = logging.getLogger(__name__)

    config = get_config("S1")

    # Read the input FASTA file
    aa_list = read_fasta(config["input_file"])

    # Generate fragments as AA sequences
    sorted_fragments_14aa = generate_fragments(config["wSet"], aa_list, 14, 1)
    sorted_fragments_14aa_G4S = generate_fragments(config["wSet"], aa_list, 14, 3)
    sorted_fragments_14aa_A5 = generate_fragments(config["wSet"], aa_list, 14, 3)
    sorted_fragments_22aa = generate_fragments(config["wSet"], aa_list, 22, 3)

    # Add overhangs for 14aa fragments
    sorted_fragments_14aa = add_overhangs(sorted_fragments_14aa, config["14_aa_overhangs"][0], config["14_aa_overhangs"][1])
    # Add overhangs for 14aa G4S fragments
    sorted_fragments_14aa_G4S = add_overhangs(sorted_fragments_14aa_G4S, config["14_aa_G4S_overhangs"][0], config["14_aa_G4S_overhangs"][1])
    # Add overhangs for 14aa A5 fragments
    sorted_fragments_14aa_A5 = add_overhangs(sorted_fragments_14aa_A5, config["14_aa_A5_overhangs"][0], config["14_aa_A5_overhangs"][1])
    # Add overhangs for 22aa fragments
    sorted_fragments_22aa = add_overhangs(sorted_fragments_22aa, config["22_aa_overhangs"][0], config["22_aa_overhangs"][1])

    logger.info(f"Number of 14aa fragments: {len(sorted_fragments_14aa)}")
    logger.info(f"Number of 14aa G4S fragments: {len(sorted_fragments_14aa_G4S)}")
    logger.info(f"Number of 14aa A5 fragments: {len(sorted_fragments_14aa_A5)}")
    logger.info(f"Number of 22aa fragments: {len(sorted_fragments_22aa)}")

    # Merge all separate fragment lists into one complete list
    sorted_fragments = sorted_fragments_22aa + sorted_fragments_14aa + sorted_fragments_14aa_A5 + sorted_fragments_14aa_G4S
    all_fragments = [frag['AAfragment'] for frag in sorted_fragments]
    logger.info(f"Number of unique Amino Acid fragments: {len(set(set(all_fragments)))}")
    logger.info(f"Number of written fragments: {len(all_fragments)}")

    # Write the merged fragments to a file
    with open(config["output_name"], "w") as f:
        f.write("Sequence\n")
        for frag in sorted_fragments:
            f.write(f"{frag['DNAfragment']}\n")
    logger.info(f"Written fragments to file: {config['output_name']}")
    logger.info(f"Execution time: {datetime.now() - start_time} seconds")


if __name__ == "__main__":
    main()