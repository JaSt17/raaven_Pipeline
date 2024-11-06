# Import necessary libraries
import random
from collections import defaultdict
import pandas as pd
from Bio.Seq import Seq

inverted_codon_table = {
            'E': ['GAG', 'GAA'],
            'S': ['TCG', 'AGT', 'TCT', 'AGC', 'TCA', 'TCC'],
            'F': ['TTC', 'TTT'],
            'K': ['AAA', 'AAG'],
            'H': ['CAC', 'CAT'],
            'C': ['TGT', 'TGC'],
            'Q': ['CAA', 'CAG'],
            'R': ['CGG', 'CGA', 'AGA', 'AGG', 'CGT', 'CGC'],
            'A': ['GCT', 'GCA', 'GCG', 'GCC'],
            'I': ['ATA', 'ATC', 'ATT'],
            'P': ['CCC', 'CCG', 'CCA', 'CCT'],
            'T': ['ACA', 'ACG', 'ACC', 'ACT'],
            'W': ['TGG'],
            'N': ['AAT', 'AAC'],
            'L': ['CTC', 'CTT', 'CTG', 'TTA', 'CTA', 'TTG'],
            'V': ['GTC', 'GTT', 'GTA', 'GTG'],
            'G': ['GGG', 'GGC', 'GGA', 'GGT'],
            'M': ['ATG'],
            'Y': ['TAC', 'TAT'],
            'D': ['GAC', 'GAT']
            }


def load_wSet(file_path: str) -> pd.DataFrame:
    """ 
    Reads the wSet file and returns a pandas DataFrame with the data.

    :param file_path: The path to the wSet file
    :return: A pandas DataFrame with the data
    """
    try:
        wSet = pd.read_csv(file_path, index_col=0)
    except Exception as e:
        print(f"An error occurred: {e}")
        return None
    return wSet


def gene_codon(seq: str, wSet: pd.DataFrame, organism: str = "ec", max_opt: bool = True, scale: float = 0.5, numcode: int = 1) -> str:
    """
    Optimizes the codon usage for a given organism.
    
    :param seq: DNA sequence
    :param wSet: Codon usage DataFrame
    :param organism: Organism for codon optimization
    :param max_opt: Whether to use maximum codon optimization
    :param scale: Scaling factor for randomness in non-maximizing optimization
    :param numcode: Genetic code to use for translation
    :return: Codon-optimized DNA sequence
    """
    
    # Load codon weights
    try:
        codon_usage_weights = wSet.loc[organism].to_dict()
    except KeyError:
        raise ValueError(f"The organism {organism} is not present in the wSet DataFrame")
    
    if len(seq) == 0:
        return ""
    
    # Translate the sequence to amino acids
    amino_seq = str(Seq(seq).translate(table=numcode))
    
    optimized_seq = []

    # Codon optimization
    for aa in amino_seq:
        codons = inverted_codon_table.get(aa, [])
        if max_opt:
            # Maximum optimization: select codon with highest weight
            optimal_codon = max(codons, key=lambda c: codon_usage_weights.get(c.lower(), 0))
        else:
            num_codons = len(codons)
            if num_codons == 1:
                optimal_codon = codons[0]
            else:
                # Scale-adjusted random selection
                cutoff = max(1, round(num_codons * scale))
                top_codons = sorted(codons, key=lambda c: codon_usage_weights.get(c.lower(), 0), reverse=True)[:cutoff]
                optimal_codon = random.choice(top_codons)
        
        optimized_seq.append(optimal_codon)

    return "".join(optimized_seq)



def aatodna(in_aa: str, wSet: pd.DataFrame, species: str = "hsa", opt: bool = True, max_opt=True) -> str:
    """
    Translates an amino acid sequence to a DNA sequence using codon optimization.

    :param in_aa: Amino acid sequence
    :param wSet: Codon usage DataFrame
    :param species: Organism for codon optimization (default: human)
    :param opt: If True, optimize codon usage
    :param max_opt: Whether to use maximum codon optimization
    :return: Codon-optimized DNA sequence
    """
    if not isinstance(in_aa, str):
        raise ValueError("Input sequence must be a string")
    if in_aa == "":
        raise ValueError("Input sequence cannot be empty")

    in_aa = in_aa.upper()
    
    human_codon = {
    'A': 'gcc',
    'C': 'tgc',
    'D': 'gac',
    'E': 'gag',
    'F': 'ttc',
    'G': 'ggc',
    'H': 'cac',
    'I': 'atc',
    'K': 'aag',
    'L': 'ctg',
    'M': 'atg',
    'N': 'aac',
    'P': 'ccc',
    'Q': 'cag',
    'R': 'cgg',
    'S': 'agc',
    'T': 'acc',
    'V': 'gtg',
    'W': 'tgg',
    'Y': 'tac'
    }

    # Translate amino acid sequence to DNA sequence
    dna_seq = ''.join([human_codon[aa]for aa in in_aa]).upper()
    
    if opt:
        optimized_dna_seq = gene_codon(dna_seq, wSet, organism=species, max_opt=max_opt)
        return optimized_dna_seq.upper()

    return dna_seq.upper()
