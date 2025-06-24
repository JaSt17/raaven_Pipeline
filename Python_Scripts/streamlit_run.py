import streamlit as st
import pandas as pd
import os
import sys
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode

# ---------- Setup and Input Validation ----------

# Get --input_dir from command-line
input_dir = None
if "--input_dir" in sys.argv:
    idx = sys.argv.index("--input_dir") + 1
    if idx < len(sys.argv):
        input_dir = sys.argv[idx]

if not input_dir:
    st.error("❌ You must provide --input_dir. Run with:\n`streamlit run streamlit_run.py -- --input_dir <path>`")
    st.stop()

if not os.path.isdir(input_dir):
    st.error(f"❌ `{input_dir}` is not a valid directory.")
    st.stop()

# Setup page
plasmid_run = os.path.basename(os.path.normpath(input_dir))
library_name = os.path.basename(os.path.dirname(os.path.normpath(input_dir)))
report_name = f"{library_name}-{plasmid_run}"
st.set_page_config(page_title=report_name, layout="wide")
st.title(f"Library Report: `{report_name}`")

# ---------- Section: Library Sequencing Summary ----------
with st.expander("Library Sequencing Summary", expanded=False):
    summary_file = os.path.join(input_dir, "sequencing_summary.csv")
    if os.path.exists(summary_file):
        df = pd.read_csv(summary_file)
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Sequencing Summary")
            st.dataframe(df, height=420, use_container_width=True)
        with col2:
            pie_chart = os.path.join(input_dir, "plots/barcode_pie_chart.png")
            if os.path.exists(pie_chart):
                st.subheader("Barcode Distribution")
                st.image(pie_chart, use_container_width=True)
    else:
        st.error("❌ Sequencing summary file not found.")

    barcode_reads_file = os.path.join(input_dir, "plots/barcode_reads.svg")
    fragments_reads_file = os.path.join(input_dir, "plots/fragment_reads.svg")
    if os.path.exists(barcode_reads_file) and os.path.exists(fragments_reads_file):
        st.subheader("Barcode and Fragment Reads Logos")
        st.image(barcode_reads_file, caption="Barcode Reads", use_container_width=True)
        st.image(fragments_reads_file, caption="Fragment Reads", use_container_width=True)

    unique_barcodes_logo_file = os.path.join(input_dir, "plots/unique_barcodes_logo.svg")
    unique_fragments_logo_file = os.path.join(input_dir, "plots/unique_fragments_logo.svg")
    if os.path.exists(unique_barcodes_logo_file) and os.path.exists(unique_fragments_logo_file):
        st.subheader("Unique Barcodes and Fragments Logos")
        col1, col2 = st.columns(2)
        col1.image(unique_barcodes_logo_file, caption="Unique Barcodes Logo", use_container_width=True)
        col2.image(unique_fragments_logo_file, caption="Unique Fragments Logo", use_container_width=True)

# ---------- Section: Found Sample Barcodes Summary ----------
with st.expander("Sample Barcodes Summary", expanded=False):
    report_file = os.path.join(input_dir, "found_barcode_report.csv")
    if os.path.exists(report_file):
        df = pd.read_csv(report_file)
        df = df.rename(columns={
            "Name": "Sample Name", "BC_reads": "Found Barcode Reads", "unique_BC": "Found Unique Barcodes",
            "matched_BC_reads": "Barcode Reads matched to Library", "matched_unique_BC": " Unique Barcodes matched to Library"
        })

        st.write("This section provides a summary of the barcodes found in all samples. It includes all barcodes with a mean read count per million reads greater than 1.\n"
                "The Found in Library column indicates whether the barcode was found in the library. So if we can match the barcode to a corresponding fragment in the library, it is marked as True.")

        gb = GridOptionsBuilder.from_dataframe(df)
        gb.configure_default_column(filterable=True, sortable=True, resizable=True)
        gb.configure_selection("single", use_checkbox=True)
        gridOptions = gb.build()

        barcode_summary = AgGrid(df, gridOptions=gridOptions, height=400, theme="alpine",
                                update_mode=GridUpdateMode.SELECTION_CHANGED, fit_columns_on_grid_load=True)

        selected_rows = barcode_summary.get("selected_rows", [])
        if selected_rows is not None:
            selected_sample = selected_rows["Sample Name"][0]
            if selected_sample:
                filename = os.path.join(input_dir, f"found_barcodes/found.{selected_sample}.csv")
                if os.path.exists(filename):
                    found_df = pd.read_csv(filename)
                    found_df = found_df[["BC", "Reads", "Peptide", "RNAcount", "RNAcount_per_mil_reads"]]
                    found_df = found_df.rename(columns={
                        "BC": "Barcode",
                        "Reads": "Fragment",
                        "RNAcount": "Raw RNA count",
                        "RNAcount_per_mil_reads": "RNA count per mil reads"
                    })
                    st.subheader(f"Details for Sample: {selected_sample}")
                    st.dataframe(found_df, use_container_width=True)
                else:
                    st.warning(f"No barcode file found for sample: {selected_sample}")
    else:
        st.warning("❌ `found_barcode_report.csv` not found.")

# ---------- Section: All Found Barcodes ----------
with st.expander("Top List Barcodes", expanded=False):
    found_barcodes_file = os.path.join(input_dir, "found_barcodes.csv")
    if os.path.exists(found_barcodes_file):
        df = pd.read_csv(found_barcodes_file)
        num_samples = len(df.columns) - 5
        df = df[df['Count_per_mil_reads_mean'] > 1]
        df.rename(columns={
            "BC": "Found Barcode",
            "Found_in_Lib": "Found in Library",
            "Count_per_mil_reads_mean": "Mean Reads per Million Reads",
            "Max_per_mil_reads": "Max Reads per Million Reads in a Sample",
            "Samples_Found_In": f"Samples Found in ... out of {num_samples}"
        }, inplace=True)

        st.write("This section summarizes barcodes found across all samples.")

        gb = GridOptionsBuilder.from_dataframe(df)
        gb.configure_default_column(filterable=True, sortable=True, resizable=True)
        AgGrid(df, gridOptions=gb.build(), height=400, theme="alpine", allow_unsafe_jscode=True)
    else:
        st.warning("❌ `found_barcodes.csv` not found.")

# ---------- Section: Final Fragment Summary ----------
with st.expander("Summary of Found Fragments", expanded=False):
    st.write("""
    This section provides a summary of all fragments that were detected by matching barcodes found in the samples to the library barcodes.\n
    You can take a closer look at the fragments by selecting them from the table below.\n
    If *1* Fragment is selected, you will see the details of the fragment across all groups.\n
    If *2* Fragments are selected, you will see a comparison of the two fragments across all groups.\n
    """)
    summary_file = os.path.join(input_dir, "final_fragments_summary.csv")
    if os.path.exists(summary_file):
        df = pd.read_csv(summary_file)
        df = df[["Group", "Peptide", "LLinker", "Sequence", "RLinker", "RNAcount", "RNAcount_ratio", "BC_count", "BC", "in_subsets"]]
        df.rename(columns={
            "BC": "Barcode(s)",
            "RNAcount": "Reads per Million Reads",
            "RNAcount_ratio": "Abundance in Group",
            "BC_count": "Barcode Count",
            "in_subsets": "Fraction of groups represented in subsets"
        }, inplace=True)

        options = [g for g in df["Group"].unique() if g != "Plasmid_Library"]
        selected_groups = st.multiselect("Select Group(s):", options, default=options)

        filtered_df = df[df["Group"].isin(selected_groups)]

        gb = GridOptionsBuilder.from_dataframe(filtered_df)
        gb.configure_default_column(filterable=True, sortable=True, resizable=True)
        gb.configure_selection("multiple", use_checkbox=True)
        gridOptions = gb.build()

        final_fragments = AgGrid(filtered_df, gridOptions=gridOptions, height=600, theme="alpine",
                                update_mode=GridUpdateMode.SELECTION_CHANGED, allow_unsafe_jscode=True, fit_columns_on_grid_load=True)
        selected_rows = final_fragments.get("selected_rows", [])

        if selected_rows is not None:
            if len(selected_rows) > 2:
                st.warning("⚠️ Only the first 2 rows will be used.")
                selected_rows = selected_rows[:2]

            if len(selected_rows) == 1:
                peptide = selected_rows['Peptide'][0]
                match = df[df["Peptide"] == peptide]
                st.subheader(f"Details for Peptide: {peptide}")
                st.dataframe(match, use_container_width=True)

            elif len(selected_rows) == 2:
                p1, p2 = selected_rows['Peptide'][0], selected_rows['Peptide'][1]
                m1 = df[df["Peptide"] == p1]
                m2 = df[df["Peptide"] == p2]
                cmp = pd.merge(m1, m2, on="Group", suffixes=(f" ({p1})", f" ({p2})"))
                for col in cmp.columns:
                    if "Reads per Million Reads" in col:
                        cmp[col] = cmp[col].fillna(0)
                cmp["Ratio"] = cmp.apply(
                    lambda row: row[f"Reads per Million Reads ({p1})"] / row[f"Reads per Million Reads ({p2})"]
                    if row[f"Reads per Million Reads ({p1})"] > 0 and row[f"Reads per Million Reads ({p2})"] > 0
                    else None, axis=1
                )
                st.subheader(f"Comparison: {p1} vs {p2}")
                st.dataframe(cmp[["Group", f"Reads per Million Reads ({p1})", f"Reads per Million Reads ({p2})", "Ratio"]],
                            use_container_width=True)

# ---------- Section: Fragment Group Comparison ----------
with st.expander("Compare Fragment Abundance Across Samples/Tissues", expanded=False):
    st.write("""
        Here we can compare the effectiveness of fragments across different groups.\n
        You can select a group to compare against other groups. The first group will be shown in the first column, while the other groups will be shown in the following columns.\n
        The comparison is based on the reads per million reads for each fragment in the selected groups.\n
        """)
    df = pd.read_csv(os.path.join(input_dir, "final_fragments_summary.csv"))
    df.rename(columns={"in_subsets": "Fraction of groups represented in subsets",
                        "RNAcount": "Mean Reads per Million Reads"}, inplace=True)
    df = df[df["Fraction of groups represented in subsets"] < 1.0]
    options = df["Group"].unique()

    col1, col2 = st.columns(2)
    first_group = col1.selectbox("Group to Compare", options)
    compare_groups = col2.multiselect("Compare To:", [g for g in options if g != first_group])

    if st.button("Compare", use_container_width=True) and compare_groups:
        display_df = df[df["Group"] == first_group][["Peptide", "Mean Reads per Million Reads", "Fraction of groups represented in subsets"]]
        display_df.rename(columns={
            "Mean Reads per Million Reads": first_group,
            "Fraction of groups represented in subsets": f"Fraction in {first_group}"
        }, inplace=True)

        for group in compare_groups:
            gd = df[df["Group"] == group][["Peptide", "Mean Reads per Million Reads", "Fraction of groups represented in subsets"]]
            gd.rename(columns={
                "Mean Reads per Million Reads": group,
                "Fraction of groups represented in subsets": f"Fraction in {group}"
            }, inplace=True)
            display_df = pd.merge(display_df, gd, on="Peptide", how="left")

        display_df.sort_values(by=first_group, ascending=False, inplace=True)
        st.subheader(f"Comparison of {first_group} with {', '.join(compare_groups)}")
        st.dataframe(display_df, use_container_width=True)

# ---------- Section: Heatmaps ----------
with st.expander("Heatmap Visualization of Fragments", expanded=False):
    heatmap_dir = os.path.join(input_dir, "plots/All_variants/heatmaps/")
    if os.path.exists(heatmap_dir):
        heatmaps = [
            "Plasmid_Amino.png", "Packaged_Amino.png", "Infective_Amino.png",
            "Packaged_vs_Plasmid_Amino.png", "Infective_vs_Plasmid_Amino.png", "Infective_vs_Packaged_Amino.png"
        ]
        files = [f for f in heatmaps if os.path.exists(os.path.join(heatmap_dir, f))]
        if files:
            cols = st.columns(len(files))
            for i, f in enumerate(files):
                with cols[i]:
                    st.image(os.path.join(heatmap_dir, f), caption=f, use_container_width=True)
        else:
            st.warning("No heatmap files found.")
    else:
        st.warning("Heatmap directory not found.")
