# Import necessary libraries
import random
import pandas as pd
from Bio.Seq import Seq
import os
import re
from plotting_functions import plot_amino_acid_heatmap, plot_aa_deviation_difference, plot_quantities

# Define the inverted codon table for codon optimization
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
    
    # set the index to the organism name
    wSet.set_index('organism', inplace=True)

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


def extract_logging_info_and_write_csv(logs_folder: str, output_csv_path: str, threshold=False):
    """
    Extracts metrics from S2.log and S3.log in the specified folder and writes a summary CSV.
    
    Args:
        logs_folder (str): Path to the folder containing S2.log and S3.log.
        output_csv_path (str): Path to save the output CSV summary.
    """
    s2_path = os.path.join(logs_folder, "S2.log")
    s3_path = os.path.join(logs_folder, "S3.log")
    data = {}

    # --- Parse S2.log ---
    with open(s2_path, "r") as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if "bbduk2 barcode extraction summary:" in line:
                if "Input:" in lines[i+1]:
                    match = re.search(r'(\d+)\s+reads', lines[i+1])
                    if match:
                        data["Number of Reads"] = int(match.group(1))
                if "Result:" in lines[i+5]:
                    match = re.search(r'Result:.*\t(\d+)\s+reads\s+\(([\d.]+%)\)', lines[i+5])
                    if match:
                        data["Reads match Barcode Pattern"] = f"{match.group(1)} ({match.group(2)})"
            elif "bbduk2 fragment extraction summary:" in line:
                if "Result:" in lines[i+5]:
                    match = re.search(r'Result:.*\t(\d+)\s+reads\s+\(([\d.]+%)\)', lines[i+5])
                    if match:
                        data["Reads match Fragment Pattern"] = f"{match.group(1)} ({match.group(2)})"
            elif "Paired reads:" in line:
                match = re.search(r'Paired reads:\s+(\d+)', line)
                if match:
                    data["Paired Reads"] = int(match.group(1))

    # --- Parse S3.log ---
    with open(s3_path, "r") as f:
        lines = f.readlines()
        for line in lines:
            if "Number of unique fragment reads" in line:
                match = re.search(r'reads: (\d+)', line)
                if match:
                    data["Unique Fragments"] = int(match.group(1))
            elif "Number of unqiue barcode reads" in line:
                match = re.search(r'reads: (\d+)', line)
                if match:
                    data["Unique Barcodes"] = int(match.group(1))
            elif "Single read barcodes" in line:
                match = re.search(r'Single read barcodes\s+(\d+)\s+\(\s*([\d.]+%)\)', line)
                if match:
                    data["Single-read barcodes"] = f"{match.group(1)} ({match.group(2)})"
            elif "Definitve read barcodes" in line:
                match = re.search(r'Definitve read barcodes\s+(\d+)\s+\(\s*([\d.]+%)\)\s+(\d+)\s+\(\s*([\d.]+%)\)', line)
                if match:
                    data["Definitive barcodes (reads)"] = f"{match.group(1)} ({match.group(2)})"
                    data["Unique Definitive barcodes"] = f"{match.group(3)} ({match.group(4)})"
            elif "Chimeric_Def barcodes" in line:
                if threshold:
                    match = re.search(r'Chimeric_Def barcodes \(ratio>(\d.\d)\)\s+(\d+)\s+\(\s*([\d.]+%)\)\s+(\d+)\s+\(\s*([\d.]+%)\)', line)
                    if match:
                        data[f"Chimeric_Def Barcodes (reads) (ratio>{match.group(1)})"] = f"{match.group(2)} ({match.group(3)})"
                        data[f"Unique Chimeric_Def Barcodes (ratio>{match.group(1)})"] = f"{match.group(4)} ({match.group(5)})"
            elif "Chimeric read barcodes" in line:
                match = re.search(r'Chimeric read barcodes\s+(\d+)\s+\(\s*([\d.]+%)\)\s+(\d+)\s+\(\s*([\d.]+%)\)', line)
                if match:
                    data["Chimeric Barcodes (reads)"] = f"{match.group(1)} ({match.group(2)})"
                    data["Unique Chimeric Barcodes"] = f"{match.group(3)} ({match.group(4)})"

    # --- Save as CSV ---
    df = pd.DataFrame(list(data.items()), columns=["Metric", "Value"])
    df.to_csv(output_csv_path, index=False)
    print(f"Summary saved to {output_csv_path}")
    
    
def create_summary_plots(final_summary: pd.DataFrame, plot_dir_path: str, array_size: int):
    """
    This function combines all the plotting functions to create a summay funcion.
    It runs after all output tabels were created and saved to get an overview of the data.
    """
    # Check if the datafrmae has the coulmn "in_reference"
    if 'in_reference' in final_summary.columns:
        # If the column is present, filter the DataFrame to only include rows where 'in_reference' is True
        ref_final_summary = final_summary[final_summary['in_reference'] == True]
        prefixlist = ["All_variants", "In_reference_variants"]
        df_list = [final_summary, ref_final_summary]
    else:
        prefixlist = ["All_variants"]
        df_list = [final_summary]
        
    
    for prefix, final_summary in zip(prefixlist, df_list):
        # create the directory for the plots for each prefix
        prefix_dir = os.path.join(plot_dir_path, prefix)
        if not os.path.exists(prefix_dir):
            os.makedirs(prefix_dir)
        # crete the heatmap directory if it does not exist
        heatmap_dir = os.path.join(prefix_dir, "heatmaps")
        if not os.path.exists(heatmap_dir):
            os.makedirs(heatmap_dir)
            
        if 'Infective_AAVs' in final_summary['Group'].unique():
            df_infective = final_summary[final_summary['Group'] == "Infective_AAVs"]
            plot = plot_amino_acid_heatmap(df_infective, "Infective_AAVs", number_of_top=df_infective.shape[0])
            plot.savefig(f"{heatmap_dir}/Infective_Amino.png", dpi=600, bbox_inches='tight', transparent=True)

        if 'DNAse_resistant_AAVs' in final_summary['Group'].unique():
            # Suppress SettingWithCopyWarning by using .loc
            final_summary.loc[:, 'Group'] = final_summary['Group'].replace({'DNAse_resistant_AAVs': 'Packaged AAVs'})

            df_packaged = final_summary[final_summary['Group'] == "Packaged AAVs"]
            plot = plot_amino_acid_heatmap(df_packaged, "Packaged AAVs", number_of_top=df_packaged.shape[0])
            plot.savefig(f"{heatmap_dir}/Packaged_Amino.png", dpi=600, bbox_inches='tight', transparent=True)

            if 'Infective_AAVs' in final_summary['Group'].unique():
                plot = plot_aa_deviation_difference(final_summary, "Packaged AAVs", "Infective_AAVs")
                plot.savefig(f"{heatmap_dir}/Packaged_vs_Infective_Amino.png", dpi=600, bbox_inches='tight', transparent=True)

        if 'Plasmid_Library' in final_summary['Group'].unique():
            df_plasmid = final_summary[final_summary['Group'] == "Plasmid_Library"]
            plot = plot_amino_acid_heatmap(df_plasmid, "Plasmid_Library", number_of_top=df_plasmid.shape[0])
            plot.savefig(f"{heatmap_dir}/Plasmid_Amino.png", dpi=600, bbox_inches='tight', transparent=True)

            if 'Packaged AAVs' in final_summary['Group'].unique():
                plot = plot_aa_deviation_difference(final_summary, "Plasmid_Library", "Packaged AAVs")
                plot.savefig(f"{heatmap_dir}/Packaged_vs_Plasmid_Amino.png", dpi=600, bbox_inches='tight', transparent=True)

            if 'Infective_AAVs' in final_summary['Group'].unique():
                plot = plot_aa_deviation_difference(final_summary, "Plasmid_Library", "Infective_AAVs")
                plot.savefig(f"{heatmap_dir}/Infective_vs_Plasmid_Amino.png", dpi=600, bbox_inches='tight', transparent=True)

        
        # Create a summary quantative plot
        # check if all groups are present in the final_summary DataFrame
        if {'Plasmid_Library', 'Packaged AAVs', 'Infective_AAVs'}.issubset(final_summary['Group'].unique()):
            plot = plot_quantities(final_summary,
                        {'Plasmid_Library':'Plasmid Library', 'Packaged AAVs':'Packaged AAVs', 'Infective_AAVs': 'Infective AAVs'},
                        {"Array": array_size})
            plot.savefig(f"{prefix_dir}/quantitative_summary.png", dpi=600, bbox_inches='tight', transparent=True)