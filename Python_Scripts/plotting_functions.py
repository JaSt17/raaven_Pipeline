import os
import subprocess
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import logomaker

def generate_sequence_logo_from_fastq(input_fastq_gz_path: str, output_logo_path: str, library_name: str = ""):
    # Determine corresponding .fastq path
    input_fastq_path = input_fastq_gz_path.replace(".fastq.gz", ".fastq")

    # Step 1: Check if .fastq exists; if not, extract it
    if not os.path.isfile(input_fastq_path):
        subprocess.run(["gunzip", "-k", input_fastq_gz_path], check=True)

    base_counts = {}
    num_sequences = 0
    read_length = None  # to determine figsize later

    # Step 2: Parse sequences and count base frequencies per position
    with open(input_fastq_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)

            if read_length is None:
                read_length = len(seq)  # set read length based on first sequence

            for pos, base in enumerate(seq):
                if pos not in base_counts:
                    base_counts[pos] = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
                if base not in base_counts[pos]:
                    base_counts[pos][base] = 0
                base_counts[pos][base] += 1

            num_sequences += 1

    # Step 3: Normalize base counts to probabilities
    sorted_positions = sorted(base_counts.keys())
    bases = ['A', 'C', 'G', 'T', 'N']
    prob_matrix = []

    for pos in sorted_positions:
        total = sum(base_counts[pos].values())
        row = [base_counts[pos].get(base, 0) / total for base in bases]
        prob_matrix.append(row)

    probability_df = pd.DataFrame(prob_matrix, columns=bases)

    # Step 4: Generate logo with dynamic figsize
    if read_length is None:
        return

    logo = logomaker.Logo(probability_df, figsize=(read_length*0.5, 5))
    logo.style_glyphs(color_scheme='colorblind_safe')
    logo.style_xticks(anchor=0, spacing=1)

    plt.title(f"{library_name}", fontsize=14)
    plt.xlabel("Position in the Fragment", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)

    logo.ax.patch.set_alpha(0.0)
    logo.fig.patch.set_alpha(0.0)

    # Step 5: Save plot
    os.makedirs(os.path.dirname(output_logo_path), exist_ok=True)
    plt.savefig(output_logo_path, transparent=True)
    plt.close()

def generate_sequence_logo_from_fasta(input_fasta_path: str, output_logo_path: str, library_name: str = ""):
    base_counts = {}
    num_sequences = 0
    read_length = None  # to determine figsize later

    # Step 1: Parse sequences and count base frequencies per position
    with open(input_fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq)

            if read_length is None:
                read_length = len(seq)  # set read length based on first sequence

            for pos, base in enumerate(seq):
                if pos not in base_counts:
                    base_counts[pos] = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
                if base not in base_counts[pos]:
                    base_counts[pos][base] = 0
                base_counts[pos][base] += 1

            num_sequences += 1

    # Step 2: Normalize base counts to probabilities
    sorted_positions = sorted(base_counts.keys())
    bases = ['A', 'C', 'G', 'T', 'N']
    prob_matrix = []

    for pos in sorted_positions:
        total = sum(base_counts[pos].values())
        row = [base_counts[pos].get(base, 0) / total for base in bases]
        prob_matrix.append(row)

    probability_df = pd.DataFrame(prob_matrix, columns=bases)

    # Step 3: Generate logo with dynamic figsize
    if read_length is None:
        return

    logo = logomaker.Logo(probability_df, figsize=(read_length*0.5, 5))
    logo.style_glyphs(color_scheme='colorblind_safe')
    logo.style_xticks(anchor=0, spacing=1)

    plt.title(f"{library_name}", fontsize=14)
    plt.xlabel("Position in the Fragment", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)

    # Step 4: Save plot
    os.makedirs(os.path.dirname(output_logo_path), exist_ok=True)
    plt.savefig(output_logo_path, transparent=True)
    plt.close()