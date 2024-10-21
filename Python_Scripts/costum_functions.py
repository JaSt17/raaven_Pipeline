# Import necessary libraries
import random
import pandas as pd

human_codon = {
    'A': ['gct', 'gcc', 'gca', 'gcg'],      # Alanine
    'C': ['tgt', 'tgc'],                    # Cysteine
    'D': ['gat', 'gac'],                    # Aspartic Acid
    'E': ['gaa', 'gag'],                    # Glutamic Acid
    'F': ['ttt', 'ttc'],                    # Phenylalanine
    'G': ['ggt', 'ggc', 'gga', 'ggg'],      # Glycine
    'H': ['cat', 'cac'],                    # Histidine
    'I': ['att', 'atc', 'ata'],             # Isoleucine
    'K': ['aaa', 'aag'],                    # Lysine
    'L': ['tta', 'ttg', 'ctt', 'ctc', 'cta', 'ctg'],  # Leucine
    'M': ['atg'],                          # Methionine (start codon)
    'N': ['aat', 'aac'],                    # Asparagine
    'P': ['cct', 'ccc', 'cca', 'ccg'],      # Proline
    'Q': ['caa', 'cag'],                    # Glutamine
    'R': ['cgt', 'cgc', 'cga', 'cgg', 'aga', 'agg'],  # Arginine
    'S': ['tct', 'tcc', 'tca', 'tcg', 'agt', 'agc'],  # Serine
    'T': ['act', 'acc', 'aca', 'acg'],      # Threonine
    'V': ['gtt', 'gtc', 'gta', 'gtg'],      # Valine
    'W': ['tgg'],                          # Tryptophan
    'Y': ['tat', 'tac'],                    # Tyrosine
    '*': ['taa', 'tag', 'tga']              # Stop codons
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


def gene_codon(aa_seq: str, wSet: pd.DataFrame, organism: str = "hsa", max_opt: bool = True) -> str:
    """
    Optimizes the codon usage for a given organism.

    :param aa_seq: Amino acid sequence
    :param wSet: Codon usage DataFrame
    :param organism: Organism for codon optimization (default: human)
    :param max_opt: Whether to use maximum codon optimization
    :return: Codon-optimized DNA sequence
    """
    try:
        codon_usage_weights = wSet.loc[organism].to_dict()
    except KeyError:
        raise ValueError(f"The organism {organism} is not present in the wSet DataFrame")

    optimized_seq = []

    for aa in aa_seq:
        codons = human_codon[aa]
        if max_opt:
            optimal_codon = max(codons, key=lambda c: codon_usage_weights.get(c, 0))
        else:
            weights = [codon_usage_weights.get(c, 0) for c in codons]
            optimal_codon = random.choices(codons, weights=weights)[0]
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

    if opt:
        optimized_dna_seq = gene_codon(in_aa, wSet, organism=species, max_opt=max_opt)
        return optimized_dna_seq.upper()

    # Translate amino acid sequence to DNA sequence
    dna_seq = ''.join([human_codon[aa][0] for aa in in_aa])
    return dna_seq.upper()
