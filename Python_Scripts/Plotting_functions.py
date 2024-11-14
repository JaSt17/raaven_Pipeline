import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_rna_counts(df: pd.DataFrame, group1: str, group2: str, gene_name: str, normalize: bool = True):
    """
    Plots the Normalized_RNAcount for two groups across the relative length of a gene.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing the data.
        group1 (str): The name of the first group.
        group2 (str): The name of the second group.
        gene_name (str): The name of the gene to analyze.
        normalize (bool): Whether to normalize the data by taking the log2 of the columns. Default is True.

    Returns:
        retruns the plot
    """
    # Filter the DataFrame for the specified gene and groups
    filtered_df = df[
        (df['GeneName'] == gene_name) &
        (df['Group'].isin([group1, group2]))
    ].copy()  # Use .copy() to avoid the SettingWithCopyWarning

    if filtered_df.empty:
        print(f"No data found for GeneName '{gene_name}' with groups '{group1}' and '{group2}'.")
        return

    # Bin the AA_rel_pos into 1% increments
    bins = np.linspace(0, 1, 101)  # 0%, 1%, ..., 100%
    labels = bins[:-1] + 0.005  # Label bins at the center
    filtered_df['AA_rel_bin'] = pd.cut(filtered_df['AA_rel_pos'], bins=bins, labels=labels, include_lowest=True)

    # Aggregate the Normalized_RNAcount for each bin and group
    aggregated = filtered_df.groupby(['Group', 'AA_rel_bin'], observed=False)['Normalized_RNAcount'].sum().reset_index()
    # Pivot the data to have groups as columns
    pivot_df = aggregated.pivot(index='AA_rel_bin', columns='Group', values='Normalized_RNAcount').fillna(0)
    
    if normalize:
        # Normalize the data by taking the log2 of the columns
        pivot_df = pivot_df.apply(lambda x: np.log2(x + 1))

    # Ensure both groups have data in all bins
    pivot_df = pivot_df.reindex(labels, fill_value=0)

    # Prepare data for plotting
    x = pivot_df.index.astype(float)
    y1 = pivot_df[group1].values
    y2 = pivot_df[group2].values

    # Create the back-to-back histogram
    fig, ax = plt.subplots(figsize=(12, 5))

    # Plot group1 and group2 with specified colors and small gaps
    ax.bar(x, y1, width=0.008, align='edge', label=group1, color='#89ABE3')    # Light blue color for group1
    ax.bar(x, -y2, width=0.008, align='edge', label=group2, color='darkblue')   # Dark red color for group2

    # Customize the plot
    ax.set_xlabel('Relative Position in Gene (%)')
    if normalize:
        ax.set_ylabel('log2 of adjusted RNA Count')
    else:
        ax.set_ylabel('adjusted RNA Count')
    ax.set_title(f'{gene_name}')
    
    # Add a legend
    legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title_fontsize='large', handlelength=0, handletextpad=0)

    # Rotate and color the legend text
    legend_colors = ['#89ABE3', 'darkblue']  # Colors corresponding to the groups
    for text, color in zip(legend.get_texts(), legend_colors):
        text.set_rotation(90)
        text.set_color(color)

    ax.axhline(0, color='black', linewidth=0.5)  # Add a horizontal line at y=0
    
    # Set x-axis ticks to show percentage and range from 0 to 100%
    ax.set_xticks(np.linspace(0, 1, 11))
    ax.set_xticklabels([f'{int(tick*100)}%' for tick in ax.get_xticks()])
    ax.set_xlim(0, 1.01)

    # Adjust y-axis to be symmetric
    max_y = max(y1.max(), y2.max())
    ax.set_ylim(-max_y * 1.1, max_y * 1.1)
    
    # change the label of the y-axis to be positive
    ax.set_yticklabels([f'{abs(tick):.0f}' for tick in ax.get_yticks()])

    # Remove the grid
    ax.grid(False)

    # Remove the top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Return the plot
    return plt.show()


def plot_amino_acid_heatmap(df, group_name, structure_name):
    # Filter the DataFrame for the specified group
    group_df = df[df['Group'] == group_name]
    group_df = group_df[group_df['structure'] == structure_name]
    
    # Get only unique peptides
    group_df = group_df.drop_duplicates(subset=['Peptide'])

    # Extract the 'Peptide' sequences
    peptides = group_df['Peptide'].dropna().tolist()

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

    # Create the heatmap with a smaller and thinner color bar
    plt.figure(figsize=(12, 8))
    sns.heatmap(
        percentage_matrix,
        annot=True,
        fmt=".1f",
        cmap="coolwarm",  # Changed color scheme to 'coolwarm'
        cbar_kws={'label': 'Percentage (%)', 'shrink': 1.0, 'aspect': 15}  # Adjusted color bar size and thickness
    )
    plt.title(f'Amino Acid Usage in: {group_name}')
    plt.xlabel('Peptide Position')
    plt.ylabel('Amino Acid')
    plt.yticks(rotation=0)
    
    #return the plot
    return plt.show()

