# Import necessary libraries
import random
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
            'D': ['GAC', 'GAT'],
            'X': ['TAA', 'TAG', 'TGA']
            }


def load_wSet() -> pd.DataFrame:
    """ 
    Creates the wSet DataFrame and returns it.
    """
    # Define header (64 codons)
    header = ['organism', "aaa","aac","aag","aat","aca","acc","acg","act",
                "aga","agc","agg","agt","ata","atc","atg","att","caa","cac","cag","cat",
                "cca","ccc","ccg","cct","cga","cgc","cgg","cgt","cta","ctc","ctg","ctt",
                "gaa","gac","gag","gat","gca","gcc","gcg","gct","gga","ggc","ggg","ggt",
                "gta","gtc","gtg","gtt","taa","tac","tag","tat","tca","tcc","tcg","tct",
                "tga","tgc","tgg","tgt","tta","ttc","ttg","ttt"]

    # Define hsa row with corresponding values
    hsa = ['hsa', 0.42,0.54,0.58,0.46,0.28,0.36,0.12,0.24,0.2,0.24,0.2,0.15,0.16,
            0.48,1,0.36,0.25,0.59,0.75,0.41,0.27,0.33,0.11,0.28,0.11,0.19,0.21,
            0.08,0.07,0.2,0.41,0.13,0.42,0.54,0.58,0.46,0.23,0.4,0.11,0.26,0.25,
            0.34,0.25,0.16,0.11,0.24,0.47,0.18,0.28,0.57,0.2,0.43,0.15,0.22,0.06,
            0.18,0.52,0.55,1,0.45,0.07,0.55,0.13,0.45]

    # Create DataFrame and add hsa row
    wSet = pd.DataFrame(columns=header)
    wSet.loc[0] = hsa

    return wSet


def gene_codon(seq: str, wSet: pd.DataFrame, organism: str = "hsa", max_opt: bool = True, scale: float = 0.5, numcode: int = 1) -> str:
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
    # change the stop codon to X
    amino_seq = amino_seq.replace("*", "X")
    
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
    'Y': 'tac',
    'X': 'tga'
    }

    # Translate amino acid sequence to DNA sequence
    dna_seq = ''.join([human_codon[aa]for aa in in_aa]).upper()
    
    if opt:
        optimized_dna_seq = gene_codon(dna_seq, wSet, organism=species, max_opt=max_opt)
        return optimized_dna_seq.upper()

    return dna_seq.upper()

