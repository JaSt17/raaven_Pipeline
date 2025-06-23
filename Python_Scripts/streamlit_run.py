import streamlit as st
import pandas as pd
import plotly.express as px
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import os
import sys
import re

# --- Parse command-line arguments ---
input_dir = None
if "--input_dir" in sys.argv:
    idx = sys.argv.index("--input_dir") + 1
    if idx < len(sys.argv):
        input_dir = sys.argv[idx]

# --- Check if input_dir was provided ---
if not input_dir:
    st.error("❌ You must provide an --input_dir argument when launching the app.\n"
                "Please run the app with:\n`streamlit run streamlit_run.py -- --input_dir <path_to_input_dir>`.")
    st.stop()

# --- Validate the directory ---
if not os.path.isdir(input_dir):
    st.error(f"❌ The specified input_dir does not exist: `{input_dir}`")
    st.stop()
    
# get plasmid run name from the input directory
plasmid_run = os.path.basename(os.path.normpath(input_dir))
# get the library name from the input directory
library_name= os.path.basename(os.path.dirname(os.path.normpath(input_dir)))
report_name = f"{library_name}-{plasmid_run}"
st.set_page_config(page_title=f"{report_name}", layout="wide")
st.title(f"Report for `{report_name}`")
    
# check if all requred files are present in the input directory
required_files = [
    "final_fragments_summary.csv",
    "found_barcodes.csv",
    "sequencing_summary.csv",
    "plots/barcode_pie_chart.png",
    "plots/fragment_reads.svg",
    "plots/unique_barcodes_logo.svg",
    "plots/barcode_reads.svg",
    "plots/unique_fragments_logo.svg"]
existing_files = [f for f in required_files if os.path.exists(os.path.join(input_dir, f))]
if not existing_files:
    st.error(f"❌ The input directory `{input_dir}` does not contain any required files. "
                "Please ensure that at least one following files are present:\n"
                f"{', '.join(required_files)}")
    st.stop()
    
    
# Library Sequencing Summary section
st.header("Library Sequencing Summary")
#check if the sequencing summary file exists
sequencing_summary_file = os.path.join(input_dir, "sequencing_summary.csv")
if os.path.exists(sequencing_summary_file):
    df = pd.read_csv(sequencing_summary_file)

    # Create two columns
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Sequencing Summary")
        st.dataframe(df, height=420, use_container_width=True)
    # Load and show the image
    summary_plot_file = os.path.join(input_dir, "plots/barcode_pie_chart.png")
    if os.path.exists(summary_plot_file):
        with col2:
            st.subheader("Barcode Distribution")
            st.image(summary_plot_file, use_container_width=True)
else:
    st.write("❌ Sequencing summary file not found. Please ensure the file exists in the input directory.")
# check if the found barcodes_reads and fragments_reads files exist
barcode_reads_file = os.path.join(input_dir, "plots/barcode_reads.svg")
fragments_reads_file = os.path.join(input_dir, "plots/fragment_reads.svg")
if os.path.exists(barcode_reads_file) and os.path.exists(fragments_reads_file):
    st.subheader("Barcode and Fragment Reads Logos")
    st.image(barcode_reads_file, caption="Barcode Reads", use_container_width=True)
    st.image(fragments_reads_file, caption="Fragment Reads", use_container_width=True)
else:
    st.write("❌ Barcode reads or fragment reads files not found. Please ensure the files exist in the input directory.")
    
unique_barcodes_logo_file = os.path.join(input_dir, "plots/unique_barcodes_logo.svg")
unique_fragments_logo_file = os.path.join(input_dir, "plots/unique_fragments_logo.svg")
if os.path.exists(unique_barcodes_logo_file) and os.path.exists(unique_fragments_logo_file):
    st.subheader("Unique Barcodes and Fragments Logos")
    col1, col2 = st.columns(2)
    with col1:
        st.image(unique_barcodes_logo_file, caption="Unique Barcodes Logo", use_container_width=True)
    with col2:
        st.image(unique_fragments_logo_file, caption="Unique Fragments Logo", use_container_width=True)
else:
    st.write("❌ Unique barcodes or fragments logo files not found. Please ensure the files exist in the input directory.")

st.header("Found Sample Barcodes Summary")
# check if the found_barcodes.csv file exists
found_barcode_report_file = os.path.join(input_dir, "found_barcode_report.csv")
if os.path.exists(found_barcode_report_file):
    df = pd.read_csv(found_barcode_report_file)
    df = df.rename(columns={"Name": "Sample Name", "BC_reads": "Found Barcode Reads", "unique_BC": "Found Unique Barcodes",
                            "matched_BC_reads": "Barcode Reads matched to Library", "matched_unique_BC": " Unique Barcodes matched to Library"})
    st.subheader("Report of found Barcodes in each Sample")
    st.write("This section provides a summary of the potential barcodes found in each sample.")
    st.dataframe(df, height=420, use_container_width=True)
else:
    st.write("❌ Found barcodes file not found. Please ensure the file exists in the input directory.")
    
# check if the found_barcodes.csv file exists
found_barcodes_file = os.path.join(input_dir, "found_barcodes.csv")
if os.path.exists(found_barcodes_file):
    df = pd.read_csv(found_barcodes_file)
    num_samples = len(df.columns)-5
    df = df[df['Count_per_mil_reads_mean'] > 1]  # Filter out rows with zero reads
    df.rename(columns={"BC": "Found Barcode", "Found_in_Lib": "Found in Library", "Count_per_mil_reads_mean": "Mean Reads per Million Reads",
                        "Max_per_mil_reads": "Max Reads per Million Reads in a Sample", f"Samples_Found_In": f"Samples Found in out of ({num_samples})"}, inplace=True)
    st.subheader("Information about all found Barcodes from the Samples")
    st.write("This section provides a summary of the barcodes found in all samples. It includes all barcodes with a mean read count per million reads greater than 1.\n"
                "The `Found in Library` column indicates whether the barcode was found in the library. So if we can match the barcode to a corresponding fragment in the library, it is marked as `True`.")
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_default_column(filterable=True, sortable=True, resizable=True)
    gridOptions = gb.build()
    AgGrid(df, gridOptions=gridOptions, height=400, theme="alpine")
else:
    st.write("❌ Found barcodes file not found. Please ensure the file exists in the input directory.")

# check if the final_fragments_summary.csv file exists
final_fragments_summary_file = os.path.join(input_dir, "final_fragments_summary.csv")
if os.path.exists(final_fragments_summary_file):
    df = pd.read_csv(final_fragments_summary_file)
    # --- Prepare Data ---
    df = df[["Group", "Peptide", "LLinker", "Sequence", "RLinker", "RNAcount", "RNAcount_ratio", "BC_count", "BC", "in_subsets"]]
    df.rename(columns={
        "BC": "Barcode(s)",
        "RNAcount": "Reads per Million Reads",
        "RNAcount_ratio": "Abundance in Group",
        "BC_count": "Barcode Count",
        "in_subsets": "Fraction of groups represented in subsets"
    }, inplace=True)

    st.subheader("Summary of Found Fragments")
    st.write("""
    This section provides a summary of all fragments that were detected by matching barcodes found in the samples to the library barcodes.\n
    You can take a closer look at the fragments by selecting them from the table below.\n
    If *1* Frament is selected, you will see the details of the fragment acrross all groups.\n
    If *2* Fragments are selected, you will see a comparison of the two fragments across all groups.\n
    """)
        # Properly exclude the selected group
    options = [group for group in df["Group"].unique() if group != "Plasmid_Library"]
    
    selected_groups = st.multiselect(
            "Select Group(s):",
            options=options,
            default=options,
            key="group_selection",
            help="Select one or more groups which should be displayed in the table below."
        )

    # --- Apply Filters ---
    filtered_df = df[df["Group"].isin(selected_groups)]
    
    # --- Initialize session state for selected peptide ---
    if "selected_peptide" not in st.session_state:
        st.session_state.selected_peptide = None

    # --- Configure AgGrid ---
    gb = GridOptionsBuilder.from_dataframe(filtered_df)
    gb.configure_default_column(filterable=True, sortable=True, resizable=True)
    # Enable row selection
    gb.configure_selection(selection_mode="multiple", use_checkbox=True)
    gridOptions = gb.build()

    # --- Render Grid and Capture Selection ---
    final_fragments = AgGrid(
        filtered_df,
        gridOptions=gridOptions,
        height=700,
        theme="alpine",
        update_mode=GridUpdateMode.SELECTION_CHANGED,
        allow_unsafe_jscode=True
    )
    
    # --- Handle selection logic ---
    selected_rows = final_fragments.get("selected_rows", [])
    if final_fragments["selected_rows"] is not None:
        if len(selected_rows) > 2:
            st.warning("⚠️ You can only select up to 2 rows. Only the first 2 will be used.")
            selected_rows = selected_rows[:2]
        # if only one row is selected, update the session state with the selected peptide
        if len(selected_rows) == 1:
            new_peptide = final_fragments['selected_rows']['Peptide'][0]
            if st.session_state.selected_peptide != new_peptide:
                st.session_state.selected_peptide = new_peptide
                matching_df = df[df["Peptide"] == st.session_state.selected_peptide]
                with st.expander(f"🔍 Selected Peptide Details: {st.session_state.selected_peptide}", expanded=True):
                    st.dataframe(matching_df, use_container_width=True)
        
        if len(selected_rows) == 2:
            peptide1 = selected_rows.iloc[0]['Peptide']
            peptide2 = selected_rows.iloc[1]['Peptide']
            matching_df1 = df[df["Peptide"] == peptide1]
            matching_df2 = df[df["Peptide"] == peptide2]
            comparison_df = pd.merge(
                matching_df1,
                matching_df2,
                on="Group",
                how="outer",  # include all groups from both DataFrames
                suffixes=(f" ({peptide1})", f" ({peptide2})")
            )
            # Replace NaNs in Reads per Million Reads columns with 0
            for col in comparison_df.columns:
                if "Reads per Million Reads" in col:
                    comparison_df[col] = comparison_df[col].fillna(0)

            # keep only relevant columns for comparison
            comparison_df = comparison_df[["Group",
                                            f"Reads per Million Reads ({peptide1})", f"Reads per Million Reads ({peptide2})"]]
            # create a new column which is the ratio of the two reads per million reads columns only clac columns where both vlaues are not 0 otherwise set to NaN
            comparison_df["Ratio"] = comparison_df.apply(
                lambda row: row[f"Reads per Million Reads ({peptide1})"] / row[f"Reads per Million Reads ({peptide2})"]
                if row[f"Reads per Million Reads ({peptide1})"] > 0 and row[f"Reads per Million Reads ({peptide2})"] > 0
                else None, axis=1)
            with st.expander(f"🔍 Comparison of Peptides: {peptide1} vs {peptide2}", expanded=True):
                st.dataframe(comparison_df, use_container_width=True)
    
    st.subheader("Comparison of Fragments Across Groups")
    st.write("""
    Here we can compare the effectiveness of fragments across different groups.\n
    You can select a group to compare against other groups. The first group will be shown in the first column, while the other groups will be shown in the following columns.\n
    The comparison is based on the reads per million reads for each fragment in the selected groups.\n
    """)
    
    # Split the page into two columns
    col1, col2 = st.columns(2)

    with col1:
        # i want all groups to be selecteable which ahve atleas 1 value < 1.0 in the 'in_subsets' column
        df = df[df['Fraction of groups represented in subsets'] < 1.0]
        options = df["Group"].unique()
        first_group = st.selectbox(
            "Select a Group to compare:",
            options=options,
            index=0,
            help="Select a group which you want to compare to the performance of other groups "
        )

    # Properly exclude the selected group
    options = [group for group in df["Group"].unique() if group != first_group]

    with col2:
        compare_groups = st.multiselect(
            "Select Group(s):",
            options=options,
            default=None,
            key="compare_group_selection",
            help="Select one or more groups to compare to the chosen group."
        )

    compare_button = st.button("Compare", use_container_width=True)
    
    # if compare button is clicked, show the comparison plot
    if compare_button and first_group and compare_groups:
        # Filter the DataFrame for the selected groups
        display_fragments = df[df["Group"]== first_group]
        # keep only the columns that are relevant for the comparison
        display_fragments = display_fragments[["Peptide", "Reads per Million Reads", "Fraction of groups represented in subsets"]]
        # renmae the Reads per Million Reads column to the first group name
        display_fragments = display_fragments.rename(columns={"Reads per Million Reads": first_group,
                                                            "Fraction of groups represented in subsets": f"Fraction in {first_group}"})
        # iterate over the compare groups and add the data to the display_fragments DataFrame
        for group in compare_groups:
            group_data = df[df["Group"] == group][["Peptide", "Reads per Million Reads", "Fraction of groups represented in subsets"]]
            group_data = group_data.rename(columns={"Reads per Million Reads": group,
                                                    "Fraction of groups represented in subsets": f"Fraction in {group}"})
            display_fragments = pd.merge(display_fragments, group_data, on="Peptide", how="left")
        
        # sort the DataFrame by the first group
        display_fragments = display_fragments.sort_values(by=first_group, ascending=False)
            
        # display the comparison DataFrame
        st.subheader(f"Comparison of {first_group} with {', '.join(compare_groups)}")
        st.write("This section provides a comparison of the reads per million reads for the selected groups.\n"
                    "The first group is shown in the first column, while the other groups are shown in the following columns.")
        st.dataframe(display_fragments, use_container_width=True)