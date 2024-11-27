import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_top_counts(df: pd.DataFrame, top_n: int, group_col: str=None, y_axis: str="RNAcount"):
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
    # if more than 50 fragments, don't show the x-axis labels
    if top_n > 50:
        plt.xticks([])
    else:
        plt.xticks(rotation=45, ha='right', fontsize=6)
    plt.title(f'Top {top_n} fragments in {group_col}')
    
    return plt.show()


def plot_rna_counts(df: pd.DataFrame, group1: str, group2: str, gene_name: str, column: str, normalize: bool = False):
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


def plot_quantities(df: pd.DataFrame, groups: list, names: list, max_value: dict, step_size: int = 10000):
    """
    Creates a polar plot with concentric rings representing different groups.
    Each ring displays an arc proportional to the group's size relative to the maximum value.
    Additionally, each ring has its own set of angular grid lines and labels indicating step values.

    Parameters:
    df : pd.DataFrame
        The input DataFrame containing at least 'Group' and 'LUTnr' columns.
    groups : list
        A list of group names to include in the plot.
    names : list
        Names of the groups to be displayed in the plot.
    max_value : dict
        A dictionary where keys are group names and values are their corresponding maximum values.
        This is used to update and determine the scaling of the plot.
    step_size : int, optional
        The step size for the grid lines and labels, by default 10000.

    Returns:
    plot : matplotlib.pyplot
    """
    # Filter the dataframe to include only specified groups
    df = df.loc[df['Group'].isin(groups)]
    
    # Get unique LUTnr within each group
    df = df.groupby(['Group', 'LUTnr']).first().reset_index()
    
    # Calculate the size of each group
    group_size = dict(df.groupby('Group').size())
    
    # Change the group names to the names provided
    group_size = {names[i]: group_size[group] for i, group in enumerate(group_size)}
    
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
    cmap = plt.get_cmap("bone")
    colors = [cmap(i / (len(group_size) + 1)) for i in range(len(group_size) + 1)]
    
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
            color="white",
            rotation=text_rotation,
            rotation_mode='anchor'
        )
        
        # Increment the radius for the next ring
        radius += size  
        
    # create a list of values for the grid lines and labels
    values = [i for i in range(0, max_val, step_size)]
    
    # sort the groups in group items assending based on their values
    group_items = sorted(group_size.items(), key=lambda x: x[1])

    # create a list of values for the grid lines and labels
    radius = 1 
    for i, (_,length) in enumerate(group_items):
        # Generate values for grid lines
        temp_values = [value for value in values if value <= length]
        
        # delete all values in temp_values from values
        for value in temp_values:
            values.remove(value)
        
        # Calculate angles corresponding to these values
        angles = [(value / max_val) * 2 * np.pi for value in temp_values]
        
        # Draw grid lines and labels
        for angle, value in zip(angles, temp_values):
            # Draw the grid line from inner to outer radius of the ring
            ax.plot([angle, angle], [radius + size, radius + size + 0.05], color='grey', linewidth=0.5, linestyle='--')
            # Place the label slightly outside the outer radius
            label_radius = radius + size + 0.1 
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
                f"{value}",
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
    return plt.show()


def plot_amino_acid_heatmap(df, group_name:str=None, structure_name:str=None):
    """
    Create a heatmap showing the percentage of each amino acid at each position in the peptide sequences.

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