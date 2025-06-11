import os
import subprocess
from Bio import SeqIO
import logomaker
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

def generate_sequence_logo_from_fastq(input_fastq_gz_path: str, output_logo_path: str, library_name: str = ""):
    """
    Generates a sequence logo from a gzipped FASTQ file and saves it as an image.

    This function processes a `.fastq.gz` file by first extracting it (if necessary),
    parsing the sequences, computing base frequencies per position, normalizing these
    to probabilities, and visualizing the result as a sequence logo using Logomaker.

    Parameters:
        input_fastq_gz_path (str): Path to the input `.fastq.gz` file containing sequencing reads.
        output_logo_path (str): Path where the resulting logo image will be saved (e.g., `.png` or `.svg`).
        library_name (str, optional): Optional title to display on the sequence logo plot. Defaults to an empty string.
    """
    # Determine corresponding .fastq path
    input_fastq_path = input_fastq_gz_path.replace(".fastq.gz", ".fastq")

    # Check if .fastq exists; if not, extract it
    if not os.path.isfile(input_fastq_path):
        subprocess.run(["gunzip", "-k", input_fastq_gz_path], check=True)

    base_counts = {}
    num_sequences = 0
    read_length = None  # to determine figsize later

    # Parse sequences and count base frequencies per position
    with open(input_fastq_path, "r") as handle:
        for record in SeqIO.parse(handle, "fastq"):
            seq = str(record.seq)

            if read_length is None:
                read_length = len(seq)  # set read length based on first sequence
            # Count base frequencies at each position
            for pos, base in enumerate(seq):
                if pos not in base_counts:
                    base_counts[pos] = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
                if base not in base_counts[pos]:
                    base_counts[pos][base] = 0
                base_counts[pos][base] += 1

            num_sequences += 1

    # Normalize base counts to probabilities
    sorted_positions = sorted(base_counts.keys())
    bases = ['A', 'C', 'G', 'T', 'N']
    prob_matrix = []

    for pos in sorted_positions:
        total = sum(base_counts[pos].values())
        row = [base_counts[pos].get(base, 0) / total for base in bases]
        prob_matrix.append(row)

    probability_df = pd.DataFrame(prob_matrix, columns=bases)

    # Generate logo with dynamic figsize
    if read_length is None:
        return
    # The figuresize is set to be half the read length in width and a fixed height
    logo = logomaker.Logo(probability_df, figsize=(read_length*0.5, 5))
    logo.style_glyphs(color_scheme='colorblind_safe')
    logo.style_xticks(anchor=0, spacing=1)

    plt.title(f"{library_name}", fontsize=14)
    plt.xlabel("Position in the Fragment", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    logo.ax.patch.set_alpha(0.0)
    logo.fig.patch.set_alpha(0.0)

    # Save plot
    os.makedirs(os.path.dirname(output_logo_path), exist_ok=True)
    plt.savefig(output_logo_path, transparent=True)
    plt.close()


def generate_sequence_logo_from_fasta(input_fasta_path: str, output_logo_path: str, library_name: str = ""):
    """
    Generates a sequence logo from a FASTA file and saves it as an image.

    This function reads nucleotide sequences from a FASTA file, computes base frequencies 
    (A, C, G, T, N) at each position across all sequences, normalizes these frequencies 
    to probabilities, and visualizes them as a sequence logo using Logomaker.

    Parameters:
        input_fasta_path (str): Path to the input `.fasta` file containing nucleotide sequences.
        output_logo_path (str): Path to save the resulting logo image (e.g., `.png`, `.svg`).
        library_name (str, optional): Optional title for the sequence logo plot. Defaults to an empty string.
    """
    base_counts = {}
    num_sequences = 0
    read_length = None  # to determine figsize later

    # Parse sequences and count base frequencies per position
    with open(input_fasta_path, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seq = str(record.seq)

            if read_length is None:
                read_length = len(seq)  # set read length based on first sequence
            # Count base frequencies at each position
            for pos, base in enumerate(seq):
                if pos not in base_counts:
                    base_counts[pos] = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
                if base not in base_counts[pos]:
                    base_counts[pos][base] = 0
                base_counts[pos][base] += 1

            num_sequences += 1

    # Normalize base counts to probabilities
    sorted_positions = sorted(base_counts.keys())
    bases = ['A', 'C', 'G', 'T', 'N']
    prob_matrix = []

    for pos in sorted_positions:
        total = sum(base_counts[pos].values())
        row = [base_counts[pos].get(base, 0) / total for base in bases]
        prob_matrix.append(row)

    probability_df = pd.DataFrame(prob_matrix, columns=bases)

    # Generate logo with dynamic figsize
    if read_length is None:
        return

    logo = logomaker.Logo(probability_df, figsize=(read_length*0.5, 5))
    logo.style_glyphs(color_scheme='colorblind_safe')
    logo.style_xticks(anchor=0, spacing=1)

    plt.title(f"{library_name}", fontsize=14)
    plt.xlabel("Position in the Fragment", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)

    # ave plot
    os.makedirs(os.path.dirname(output_logo_path), exist_ok=True)
    plt.savefig(output_logo_path, transparent=True)
    plt.close()
    
    
def plot_amino_acid_heatmap(df, group_name: str = None, structure_name: str = None, number_of_top: int = 100):
    """
    Create a heatmap showing the deviation of amino acid frequency from a uniform distribution at each peptide position.

    Parameters:
    df : pd.DataFrame
        The input DataFrame containing the peptide sequences.
    group_name : str, optional
        The name of the group to filter the data. Default is None.
    structure_name : str, optional
        The name of the structure to filter the data. Default is None.

    Returns:
    plot : matplotlib.pyplot
    """
    if group_name is not None:
        temp_df = df[df['Group'] == group_name]
    if structure_name is not None:
        temp_df = df[df['Structure'] == structure_name]

    # Get only unique peptides
    temp_df = temp_df.drop_duplicates(subset=['Peptide'])

    # Extract the 'Peptide' sequences
    peptides = temp_df['Peptide'].dropna().tolist()

    # Determine the maximum peptide length
    max_length = max(len(peptide) for peptide in peptides)

    # Standard amino acids
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')

    # Initialize a DataFrame to hold counts
    count_matrix = pd.DataFrame(0, index=amino_acids, columns=range(1, max_length + 1))

    # Count occurrences of each amino acid at each position
    for peptide in peptides:
        for position, amino_acid in enumerate(peptide, start=1):
            if amino_acid in amino_acids:
                count_matrix.at[amino_acid, position] += 1

    # Convert counts to percentages
    total_counts = count_matrix.sum(axis=0)
    percentage_matrix = count_matrix.divide(total_counts, axis=1) * 100

    # Compute deviation from uniform distribution (5%)
    deviation_matrix = percentage_matrix - 5.0

    # Use TwoSlopeNorm to center 0 in the colormap
    center = 0
    vmin = -5
    vmax = 10
    divnorm = TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)
    # reverse the colormap 
    cmap = plt.get_cmap('RdBu')
    cmap = cmap.reversed()
    
    # Create the heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        deviation_matrix,
        annot=False,
        cmap=cmap,
        norm=divnorm,  # Apply normalization here
        square=True,
        linewidths=0.003,
        linecolor='white',
        cbar_kws={'label': 'Deviation from Uniform (%)', 'shrink': 1.0, 'aspect': 50}
    )

    if group_name:
        group_name = group_name.replace("_", " ")
    plt.title(f'{group_name}\nN={number_of_top}', fontsize=12)
    plt.xlabel('Peptide Position')
    plt.ylabel('Amino Acid')
    plt.yticks(rotation=0)

    return plt


def plot_aa_deviation_difference(df, group_name_1: str = None, group_name_2: str = None):
    """
    Create a heatmap comparing the change in amino acid deviation from uniform distribution between two datasets.

    Parameters:
    df1, df2 : pd.DataFrame
        Two input DataFrames containing peptide sequences.
    group_name : str, optional
        Group filter to apply to both DataFrames.
    structure_name : str, optional
        Structure filter to apply to both DataFrames.
    number_of_top : int
        Number to show in the plot title for context (e.g., number of sequences).

    Returns:
    plot : matplotlib.pyplot
    """
    def compute_deviation(df, group_name):
        if group_name is not None:
            df = df[df['Group'] == group_name]

        df = df.drop_duplicates(subset=['Peptide'])
        peptides = df['Peptide'].dropna().tolist()
        max_length = max(len(p) for p in peptides)
        amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
        count_matrix = pd.DataFrame(0, index=amino_acids, columns=range(1, max_length + 1))

        for peptide in peptides:
            for position, aa in enumerate(peptide, start=1):
                if aa in amino_acids:
                    count_matrix.at[aa, position] += 1

        total_counts = count_matrix.sum(axis=0)
        percentage_matrix = count_matrix.divide(total_counts, axis=1) * 100
        deviation_matrix = percentage_matrix - 5.0
        return deviation_matrix

    # Compute deviation matrices
    dev1 = compute_deviation(df, group_name_1)
    dev2 = compute_deviation(df, group_name_2)

    # Align columns and rows (in case lengths or AAs are mismatched)
    dev1, dev2 = dev1.align(dev2, join='outer', fill_value=0)

    # Compute difference: how much did the deviation change from df1 to df2
    delta_matrix = dev2 - dev1

    # Setup color scaling: -2.5 to +5
    vmin, vmax, center = -2.5, 2.5, 0
    norm = TwoSlopeNorm(vmin=vmin, vcenter=center, vmax=vmax)

    # Reverse colormap so red = increase, blue = decrease
    cmap = plt.get_cmap("RdBu").reversed()

    # Plot the difference heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(
        delta_matrix,
        cmap=cmap,
        norm=norm,
        square=True,
        linewidths=0.003,
        linecolor='white',
        cbar_kws={'label': 'Change in Deviation (%)', 'shrink': 1.0, 'aspect': 50}
    )

    if group_name_1:
        group_name_1 = group_name_1.replace("_", " ")
    if group_name_2:
        group_name_2 = group_name_2.replace("_", " ")
    plt.title(f'Change in AA Deviation\n{group_name_1} vs. {group_name_2}', fontsize=12)
    plt.xlabel('Peptide Position')
    plt.ylabel('Amino Acid')
    plt.yticks(rotation=0)

    return plt


def plot_quantities(df: pd.DataFrame, groups: dict, max_value: dict, step_size: int = 10000):
    """
    Creates a polar plot with concentric rings representing different groups.
    Each ring displays an arc proportional to the group's size relative to the maximum value.
    Additionally, each ring has its own set of angular grid lines and labels indicating step values.

    Parameters:
    df : pd.DataFrame
        The input DataFrame containing at least 'Group' and 'LUTnr' columns.
    groups : dict
        A dictionary where keys are group names in the DataFrame and values are the names to display on the plot.
    max_value : dict
        A dictionary where keys are group names and values are their corresponding maximum values.
        This is used to update and determine the scaling of the plot.
    step_size : int, optional
        The step size for the grid lines and labels, by default 10000.

    Returns:
    plot : matplotlib.pyplot
    """
    # Filter the dataframe to include only specified groups
    df = df.loc[df['Group'].isin(groups.keys())]
    
    # Get unique LUTnr within each group
    df = df.groupby(['Group', 'LUTnr']).first().reset_index()
    
    # Calculate the size of each group
    group_size = dict(df.groupby('Group').size())
    
    # Change the group names to the names provided
    group_size = {groups[key]: value for key, value in group_size.items()}
    
    # Update the group sizes with the max_value provided
    group_size.update(max_value)
    
    # Determine the maximum value for scaling
    max_val = max(group_size.values())
    
    # Calculate the size (thickness) of each ring
    size = 1 / (len(groups) + 1)
    
    # Create a figure and axis with polar projection
    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})
    
    # Set zero at the top and angles to increase clockwise
    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    
    # Disable the default grid and ticks
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    # remove the outside border
    ax.spines['polar'].set_visible(False)
    
    # Generate colors using the 'bone' colormap
    cmap = plt.get_cmap("mako")
    colors = [cmap(0.1 + (i / (len(group_size) + 1)*0.8)) for i in range(len(group_size) + 1)]
    
    # Sort group sizes for consistent plotting (largest first)
    group_items = sorted(group_size.items(), key=lambda x: x[1], reverse=True)
    
    # Initialize the radius for the first ring
    radius = size
    for i, (label, length) in enumerate(group_items):
        # Calculate the angle filled
        angle_filled = (length / max_val) * 2 * np.pi  # Angle in radians
        
        # Draw the arc using ax.bar()
        ax.bar(
            x=0,                # Start angle (0 radians)
            height=size,        # Thickness of the ring
            width=angle_filled, # Angle spanned by the arc
            bottom=radius,      # Starting radius
            color=colors[i],
            edgecolor='w',
            align='edge',
        )
        # Get the position for annotation
        if angle_filled / (2*np.pi) > 0.6:
            theta = np.pi
        else:
            theta = angle_filled / 2  # Middle of the filled arc
        x = radius + size / 2     # Mid-radius of the ring
        
        # if the label starts with a * than set color to black
        if label.startswith('*'):
            # remove the * from the label
            label = label[1:]
            color = 'black'
            theta = theta + 0.15
        else:
            color = 'white'
        
        # Determine text rotation for better readability
        theta_deg = np.degrees(theta)
        if theta_deg > 180:
            text_rotation = theta_deg - 180
        else:
            text_rotation = 180 - theta_deg
        text_rotation = text_rotation if text_rotation < 90 else text_rotation - 180
        
        # Annotate with group name and size
        ax.text(
            theta,
            x,
            f"{label}\n{length}",
            ha="center",
            va="center",
            fontsize=10,
            color=color,
            rotation=text_rotation,
            rotation_mode='anchor'
        )
        
        # Increment the radius for the next ring
        radius += size  
    
    values = [max_val*i/10 for i in range(0,10,1)]
    
    # sort the groups in group items assending based on their values
    group_items = sorted(group_size.items(), key=lambda x: x[1])

    # create a list of values for the grid lines and labels
    radius = 1 
    for i, (_,length) in enumerate(group_items):
        # Generate values for grid lines
        temp_values = [value for value in values if value <= (length+max_val*0.004)]
        
        # delete all values in temp_values from values
        for value in temp_values:
            values.remove(value)
        
        # Calculate angles corresponding to these values
        angles = [(value / max_val) * 2 * np.pi for value in temp_values]
        
        # Draw grid lines and labels
        for angle, value in zip(angles, temp_values):
            # Draw the grid line from inner to outer radius of the ring
            ax.plot([angle, angle], [radius + size, radius + size + 0.02], color='grey', linewidth=1, linestyle='-')
            # Place the label slightly outside the outer radius
            label_radius = radius + size + 0.05
            # Calculate position for annotation so it is centered on the grid line
            angle_deg = np.degrees(angle)
            rotation = 360 - angle_deg if angle_deg > 180 else -angle_deg
            if abs(rotation) >= 270 or abs(rotation)<= 90:
                rotation = rotation
            else:
                rotation = rotation + 180
            # write the value on the grid line
            ax.text(
                angle,
                label_radius,
                f"{value/max_val*100:.0f}%",
                ha='center',
                va='center',
                fontsize=8,
                color='black',
                rotation=rotation,
                rotation_mode='anchor'
            )
        # break the loop if there are no more values
        if len(temp_values) == 0:
            break
        
        radius -= size  
        
    # add a title at the bottom of the plot
    ax.text(0, 0, 'Fragment Quantity per Group', ha='center', va='center', fontsize=12, color='black')
    
    # return the plot
    return plt