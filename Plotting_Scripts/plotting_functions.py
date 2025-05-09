import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from matplotlib.colors import TwoSlopeNorm

def plot_top_counts(df: pd.DataFrame, top_n: int, group_col: str=None, y_axis: str="RNAcount", y_label: str=None):
    """
    Plots the top N fragments based on a specified column.
    
    Parameters:
        df (pd.DataFrame): The input DataFrame containing the data.
        top_n (int): The number of top fragments to plot.
        group_col (str): The name of the group to filter the data. Default is None.
        y_axis (str): The column to use for the y-axis. Default is "RNAcount
    
    Returns:
        plot : matplotlib.pyplot
    """
    if group_col is not None:
        df = df[df['Group'] == group_col]
    
    df_top_n = df.sort_values(by=y_axis, ascending=False).head(top_n)
    
    plt.figure(figsize=(10, 5))
    if group_col is not None:
            sns.barplot(x='LUTnr', y=y_axis, data=df_top_n)
    else:
        sns.barplot(x='LUTnr', y=y_axis, hue='Group', data=df_top_n)
    # if more than 25 fragments, don't show the x-axis labels
    if top_n > 25:
        plt.xticks([])
    else:
        plt.xticks(rotation=45, ha='right', fontsize=6)
    plt.title(f'Top {top_n} fragments in {group_col}')
    
    sns.despine()
    
    if y_label is not None:
        plt.ylabel(y_label)
    
    
    return plt


def plot_venn_diagram(set1, set2, set3, labels=("Set 1", "Set 2", "Set 3"), title="Venn Diagram"):
    """
    Plots a Venn diagram for three sets.
    
    Parameters:
    - set1, set2, set3: The three sets of sequences.
    - labels: Tuple of three strings representing the labels for the sets.
    """
    # Calculate the subset sizes
    subsets = {
        '100': len(set1 - set2 - set3),  # Only in set1
        '010': len(set2 - set1 - set3),  # Only in set2
        '001': len(set3 - set1 - set2),  # Only in set3
        '110': len(set1 & set2 - set3),  # In set1 and set2, but not in set3
        '101': len(set1 & set3 - set2),  # In set1 and set3, but not in set2
        '011': len(set2 & set3 - set1),  # In set2 and set3, but not in set1
        '111': len(set1 & set2 & set3)   # In all three sets
    }

    # Create the Venn diagram
    plt.figure(figsize=(6, 6))
    
    # Generate colors using the colormap
    cmap = plt.get_cmap("twilight_shifted")
    colors = [cmap(0.1 + (i / (len(labels) + 1) * 0.9)) for i in range(len(labels))]
    venn = venn3(subsets=subsets, set_labels=labels, set_colors=colors)

    # Customize text properties
    for text in venn.set_labels:
        if text:
            text.set_fontsize(12)
    for text in venn.subset_labels:
        if text:
            text.set_fontsize(8)

    # Show the plot
    plt.title(title, fontsize=14)
    plt.show()
    

def plot_rna_counts(df: pd.DataFrame, group1: str, group2: str, gene_name: str, column: str, y_label: str=None, normalize: bool=False) -> plt:
    """
    Plots the Normalized_RNAcount for two groups across the relative length of a gene.

    Parameters:
        df (pd.DataFrame): The input DataFrame containing the data.
        group1 (str): The name of the first group.
        group2 (str): The name of the second group.
        gene_name (str): The name of the gene to analyze.
        column (str): The column that should be compared between the two groups.
        normalize (bool): Whether to normalize the data by taking the log2 of the columns.

    Returns:
        plot : matplotlib.pyplot
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
    aggregated = filtered_df.groupby(['Group', 'AA_rel_bin'], observed=False)[column].sum().reset_index()
    # Pivot the data to have groups as columns
    pivot_df = aggregated.pivot(index='AA_rel_bin', columns='Group', values=column).fillna(0)

    # Ensure both groups have data in all bins
    pivot_df = pivot_df.reindex(labels, fill_value=0)

    # Prepare data for plotting
    x = pivot_df.index.astype(float)
    y1 = pivot_df[group1].values
    y2 = pivot_df[group2].values
    
    # Normalize the data if specified
    if normalize:
        y1 = np.log10(y1 + 1)
        y2 = np.log10(y2 + 1)
        y_label = f'Log₁₀ {y_label}'

    # Create the back-to-back histogram
    fig, ax = plt.subplots(figsize=(12, 5))

    # Plot group1 and group2 with specified colors and small gaps
    ax.bar(x, y1, width=0.008, align='edge', label=group1, color='#89ABE3')    # Light blue color for group1
    ax.bar(x, -y2, width=0.008, align='edge', label=group2, color='darkblue')   # Dark red color for group2

    # Customize the plot
    ax.set_xlabel('Relative Position in Gene (%)')
    if y_label is not None:
        ax.set_ylabel(y_label)
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
    ax.set_xticklabels([f'{int(tick*100)}' for tick in ax.get_xticks()])
    ax.set_xlim(0, 1.01)
    # change font size of x-axis ticks
    ax.tick_params(axis='x', labelsize=8)
    ax.tick_params(axis='y', labelsize=8)

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


def create_grouped_barplot(df, tissue_col, count_col, library_col):
    """
    Creates a grouped bar plot from a DataFrame. Shows the count values of found fragments in each group, split by library.

    Parameters:
    df : pandas.DataFrame
        DataFrame containing the data about the tissue groups, count values, and libraries.
        
    tissue_col : str
        Name of the column in the DataFrame that contains the groups.
        
    count_col : str
        Name of the column in the DataFrame that contains the count values.
        
    library_col : str
        Name of the column in the DataFrame that contains the libraries.

    Returns:
    plot : matplotlib.pyplot
    """
    # Pivot the DataFrame for easier plotting
    pivot_df = df.pivot(index=tissue_col, columns=library_col, values=count_col)
    
    # Sort rows by total counts across all libraries
    pivot_df['Total'] = pivot_df.sum(axis=1)
    pivot_df = pivot_df.sort_values(by='Total', ascending=False).drop(columns='Total')
    
    # Plot setup
    libraries = pivot_df.columns
    x = np.arange(len(pivot_df))  # x locations for the groups
    width = 0.2  # Width of each bar
    spacing = 0.05  # Space between bars within a group

    # Generate colors using the colormap
    cmap = plt.get_cmap("twilight_shifted")
    colors = [cmap(0.2 + (i / (len(libraries) + 1) * 0.8)) for i in range(len(libraries))]

    fig, ax = plt.subplots(figsize=(10, 6))

    # Create bars for each library
    for i, (library, color) in enumerate(zip(libraries, colors)):
        ax.bar(x + i * (width + spacing), pivot_df[library], width, label=library, color=color)

    # Customization
    ax.set_ylabel(f'Number of Fragments (Log₁₀ Scale)', fontsize=8)
    ax.set_title('Distribution of Found Fragments Across Tissue Types and Sorted Cell Groups', fontsize=12)
    ax.set_xticks(x + (len(libraries) - 1) * (width + spacing) / 2) 
    ax.set_xticklabels(pivot_df.index, fontsize=8, rotation=45, ha='right')  # Rotate labels at 45 degrees

    ax.legend(title=library_col, fontsize=10)

    # Set y-axis to log2 scale
    ax.set_yscale('log', base=10)

    # Remove the top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    return plt


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


def plot_aa_deviation_difference(df, group_name_1: str = None, group_name_2: str = None, structure_name: str = None, number_of_top: int = 100):
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