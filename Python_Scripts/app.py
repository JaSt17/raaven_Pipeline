import streamlit as st
import pandas as pd
import os
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode
import matplotlib.pyplot as plt
import numpy as np
from functools import reduce

# ------------------ CACHE HELPERS ------------------

@st.cache_data(show_spinner=False)
def load_csv_file(path):
    return pd.read_csv(path)

@st.cache_data(show_spinner=False)
def load_annotation(annotation_file):
    df = pd.read_csv(annotation_file)
    df["Sample"] = df["Sample"].str.replace(".fastq.gz", "", regex=False)
    return df

def add_group_column(df: pd.DataFrame, input_dir: str, name_column: str = "Name") -> pd.DataFrame:
    parent_3_dir = os.path.dirname(os.path.dirname(os.path.dirname(input_dir)))
    # get the library name form the last part of the path
    library_name = os.path.basename(os.path.dirname(input_dir))
    # if the library name starts with p034, p035 or p036 set prefix to GFP
    if library_name.startswith(("p034", "p035", "p036")):
        suffix = "GFP"
    elif library_name.startswith(("p037", "p038", "p039")):
        suffix = "mCherry"
    else:
        suffix = None
    if suffix is None:
        annotation_file = os.path.join(parent_3_dir, "Seq_Data", "Samples", "annotation.csv")
    else:
        annotation_file = os.path.join(parent_3_dir, "Seq_Data", "Samples", f"annotation_{suffix}.csv")

    if not os.path.exists(annotation_file):
        raise FileNotFoundError(f"Annotation file not found at {annotation_file}")

    annotation_df = load_annotation(annotation_file)

    merged_df = df.merge(annotation_df, left_on=name_column, right_on="Sample", how="left")

    if "Group" not in merged_df.columns:
        raise ValueError("Column 'Group' not found in merged annotation data.")

    group_col = merged_df.pop("Group")
    merged_df.insert(1, "Group", group_col)

    merged_df.drop(columns=["Sample"], inplace=True)

    return merged_df

def get_file_path(input_dir, filename):
    return os.path.join(input_dir, filename)

def get_df_from_session_or_file(key, path):
    if key not in st.session_state:
        if os.path.exists(path):
            st.session_state[key] = load_csv_file(path)
        else:
            st.session_state[key] = None
    return st.session_state[key]

@st.cache_data
def load_barcode_summary(libraries_dir):
    summary_data = []

    libraries = [
        d for d in os.listdir(libraries_dir)
        if os.path.isdir(os.path.join(libraries_dir, d)) and d != "Undetermined"
    ]

    for lib in libraries:
        library_dir = os.path.join(libraries_dir, lib)
        runs = sorted([
            d for d in os.listdir(library_dir)
            if os.path.isdir(os.path.join(library_dir, d)) and d != "mRNA_barcodes"
        ])

        if not runs:
            continue

        run_data = {}
        total_num_barcodes = 0

        for run in runs:
            run_path = os.path.join(library_dir, run)
            barcode_file = get_file_path(run_path, "found_barcodes.csv")

            if not os.path.exists(barcode_file):
                continue

            try:
                df = pd.read_csv(barcode_file, usecols=["BC", "Found_in_Lib", "Count_per_mil_reads_mean"])
            except Exception as e:
                st.warning(f"Error reading {barcode_file}: {e}")
                continue

            df = df[df["Count_per_mil_reads_mean"] > 1]
            if total_num_barcodes == 0:
                total_num_barcodes = df["BC"].nunique()
            run_data[run] = df[df["Found_in_Lib"] == True]["BC"].nunique()

        summary_data.append({
            "Library": lib,
            "Total mRNA Barcodes (BC per mil > 1)": total_num_barcodes,
            **run_data
        })

    summary_df = pd.DataFrame(summary_data)
    return summary_df.sort_values(by="Library")

# ------------------ MAIN APP ------------------

PROJECTS_DIR = "Projects"
st.set_page_config(page_title="rAAVen Sequencing Analyze Tool", layout="wide")

if 'page' not in st.session_state:
    st.session_state['page'] = 'main'

def go_back():
    st.session_state.page = 'main'
    # Clear session state for report page
    keys_to_clear = ['summary_df', 'barcode_report_df', 'unique_barcodes_df', 'final_fragments_df', 'annotated_report_df', 'best_fragments_df', 'found_barcodes_df']
    for key in keys_to_clear:
        if key in st.session_state:
            del st.session_state[key]

def go_to_report(project, library, run, input_dir):
    st.session_state.project = project
    st.session_state.library = library
    st.session_state.run = run
    st.session_state.input_dir = input_dir
    st.session_state.page = 'report'

if st.session_state.page == 'main':
    st.title("rAAVen Sequencing Analyze Tool")
    st.subheader("Detailed Analysis of Library Sequencing Runs")
    st.write("This tool allows you to analyze the library sequencing data and visualize the barcode_coverage_per_samples for different projects.")

    if not os.path.exists(PROJECTS_DIR):
        st.error(f"❌ `{PROJECTS_DIR}` directory does not exist.")
        st.stop()

    projects = [d for d in os.listdir(PROJECTS_DIR) if os.path.isdir(os.path.join(PROJECTS_DIR, d))]
    col1, col2, col3 = st.columns([3, 3, 3])

    with col1:
        selected_project = st.selectbox(
            "Project",
            projects,
            index=projects.index(st.session_state.get('project', projects[0])),
            help="Select a project to analyze its library sequencing data.",
            key="main_project"
        )

    project_dir = os.path.join(PROJECTS_DIR, selected_project)
    libraries_dir = os.path.join(project_dir, "Libraries")
    libraries = [d for d in os.listdir(libraries_dir) if os.path.isdir(os.path.join(libraries_dir, d))]
    # sort the libraries alphabetically
    libraries.sort()
    for lib in ["Undetermined"]:
        if lib in libraries:
            libraries.remove(lib)
            libraries.append(lib)

    with col2:
        selected_library = st.selectbox(
            "Library",
            libraries,
            index=libraries.index(st.session_state.get('library', libraries[0])),
            help="Select a library to analyze its sequencing data.",
            key="main_library"
        )

    library_dir = os.path.join(libraries_dir, selected_library)
    runs = [d for d in os.listdir(library_dir) if os.path.isdir(os.path.join(library_dir, d))]
    # sort the runs alphabetically
    runs.sort()
    for run in ["mRNA_barcodes"]:
        if run in runs:
            runs.remove(run)
        
    with col3:
        selected_run = st.selectbox(
            "Run",
            runs,
            index=runs.index(st.session_state.get('run', runs[0])),
            help="Select a run to analyze its sequencing data.",
            key="main_run"
        )

    if st.button("Generate Report"):
        input_dir = os.path.join(library_dir, selected_run)
        go_to_report(selected_project, selected_library, selected_run, input_dir)
        st.rerun()

# ------------------ Barcode Coverage ------------------ 

    st.subheader("Combined Library Analysis")
    st.write(
        "This section allows you to combine the findings from different library runs.\n"
    )

    col1, col2, col3 = st.columns([1,1,4])
    with col1:
        selected_project_2 = st.selectbox(
            "Project",
            projects,
            index=projects.index(st.session_state.get('project', projects[0])),
            help="Select a project to analyze its library sequencing data.",
            key="barcode_project"
        )

    project_dir = os.path.join(PROJECTS_DIR, selected_project_2)
    libraries_dir = os.path.join(project_dir, "Libraries")
    libraries = [d for d in os.listdir(libraries_dir) if os.path.isdir(os.path.join(libraries_dir, d))]
    # sort the libraries alphabetically
    libraries.sort()
    for lib in ["Undetermined"]:
        if lib in libraries:
            libraries.remove(lib)
            libraries.append(lib)

    with col2:
        selected_library_2 = st.selectbox(
            "Library",
            libraries,
            index=libraries.index(st.session_state.get('library', libraries[0])),
            help="Select a library to analyze its sequencing data.",
            key="barcode_library"
        )

    library_dir = os.path.join(libraries_dir, selected_library_2)
    runs = [d for d in os.listdir(library_dir) if os.path.isdir(os.path.join(library_dir, d))]
    runs.sort()
    for run in ["mRNA_barcodes"]:
        if run in runs:
            runs.remove(run)

    # Initial default setup (only once)
    if "previous_runs" not in st.session_state:
        st.session_state.previous_runs = []
    # Display the multiselect
    with col3:
        selected_runs = st.multiselect(
            "Runs for Barcode Coverage",
            runs,
            help="Select the runs to include in the barcode coverage check.",
            key="barcode_runs"
        )

    if selected_runs != st.session_state.previous_runs:
        st.session_state.previous_runs = selected_runs  # update state

        if not selected_runs:
            st.warning("Please select at least one run to analyze.")
            st.session_state.found_barcodes_df = None
        else:
            input_dirs = [
                os.path.join(library_dir, run)
                for run in selected_runs
                if os.path.isdir(os.path.join(library_dir, run))
            ]

            found_barcodes_df_list = []
            final_fragments_df_list = []

            report_path = get_file_path(input_dirs[0], "found_barcode_report.csv")
            if not os.path.exists(report_path):
                st.warning(f"❌ `found_barcode_report.csv` not found in {input_dirs[0]}. Please run the barcode extraction step first.")
                st.session_state.found_barcodes_df_main = None
                st.session_state.final_fragments_df = None
                st.session_state.barcode_coverage_per_sample = None
                st.stop()

            report_df = pd.read_csv(report_path)

            # Efficient loading of all CSVs
            for input_dir in input_dirs:
                lib_id = os.path.basename(input_dir)

                fbc_path = get_file_path(input_dir, "found_barcodes.csv")
                frg_path = get_file_path(input_dir, "final_fragments_summary.csv")

                fbc_df = pd.read_csv(fbc_path)
                fbc_df["Found_in_Lib"] = np.where(fbc_df["Found_in_Lib"], lib_id, "Not Found")
                found_barcodes_df_list.append(fbc_df)

                frg_df = pd.read_csv(frg_path)
                frg_df["Library_ID"] = lib_id
                final_fragments_df_list.append(frg_df)

            # Combine all Found_in_Lib columns efficiently
            if found_barcodes_df_list:
                def merge_found_in_lib(df1, df2):
                    combined = df1.copy()
                    combined["Found_in_Lib"] = df1["Found_in_Lib"].combine(df2["Found_in_Lib"], lambda x, y: (
                        y if x == "Not Found" and y != "Not Found" else
                        x if y == "Not Found" else
                        [x, y] if x != y else x
                    ))
                    return combined

                found_barcodes_df = reduce(merge_found_in_lib, found_barcodes_df_list)

                # Filtering and renaming
                found_barcodes_df.sort_values(by='Count_per_mil_reads_mean', ascending=False, inplace=True)
                num_samples = len(found_barcodes_df.columns) - 5
                found_barcodes_df = found_barcodes_df[found_barcodes_df['Count_per_mil_reads_mean'] > 1]

                rename_map = {
                    "BC": "Found Barcode",
                    "Found_in_Lib": "Found in Library",
                    "Count_per_mil_reads_mean": "Mean Reads per Million Reads",
                    "Max_per_mil_reads": "Max Reads per Million Reads in a Sample",
                    "Samples_Found_In": f"Samples Found in ... out of {num_samples}"
                }
                found_barcodes_df.rename(columns=rename_map, inplace=True)
                st.session_state.found_barcodes_df_main = found_barcodes_df
                
                # create a Barcode_coverage_per_sample.csv
                sample_columns = found_barcodes_df.columns[5:]
                barcode_coverage_per_sample = pd.DataFrame({'Sample': sample_columns})
                # Total barcodes found per sample (values > 0)
                barcode_coverage_per_sample['Barcodes'] = [ (found_barcodes_df[col] > 0).sum() for col in sample_columns ]
                # Annotated barcodes found per sample (values > 0 and Found_in_Lib == True)
                barcode_coverage_per_sample['Annotated_Barcodes'] = [
                    ((found_barcodes_df[col] > 0) & (found_barcodes_df['Found in Library'])).sum()
                    for col in sample_columns]
                barcode_coverage_per_sample['Coverage in %'] = (
                    barcode_coverage_per_sample['Annotated_Barcodes'] / barcode_coverage_per_sample['Barcodes'] * 100
                ).fillna(0)
                # covert the Coverage in % to to float with 1 decimal place
                barcode_coverage_per_sample['Coverage in %'] = barcode_coverage_per_sample['Coverage in %'].round(1).astype(str) + '%'
                st.session_state.barcode_coverage_per_sample = barcode_coverage_per_sample

            if final_fragments_df_list:
                final_fragments_df = pd.concat(final_fragments_df_list, ignore_index=True)
                final_fragments_df = final_fragments_df[final_fragments_df["Group"] != "Plasmid_Library"]
                # get only unique combianations of Group, Peptide, LLinker, Sequence, RLinker and the mean of RNAcount RNAcount_ratio, append the list of BC together keep the first in_subsets
                final_fragments_df = final_fragments_df.groupby(
                    ["Group", "Peptide", "LLinker", "Sequence", "RLinker"],
                    as_index=False
                ).agg({
                    "RNAcount": "mean",
                    "RNAcount_ratio": "mean",
                    "BC_count": "sum",
                    "BC": lambda x: ', '.join(x.unique()),
                    "in_subsets": "mean",
                    "Library_ID": lambda x: ', '.join(x.unique()),
                    "Rank_in_Group": "min"
                })
                # set the BC_count to the number of unique barcodes
                final_fragments_df["BC_count"] = final_fragments_df["BC"].apply(lambda x: len(x.split(', ')))
                # set RNAcount to int
                final_fragments_df["RNAcount"] = final_fragments_df["RNAcount"].astype(int)
                # Create new ranks based on the mean RNAcount for each group
                final_fragments_df["Rank_in_Group"] = final_fragments_df.groupby("Group")["RNAcount"].rank(ascending=False, method='dense').astype(int)
                final_fragments_df.rename(columns={
                    "RNAcount": "Reads per Million Reads",
                    "RNAcount_ratio": "Abundance in Group",
                    "BC_count": "Barcode Count",
                    "BC": "Barcode(s)",
                    "in_subsets": "Fraction of groups represented in subsets",
                    "Library_ID": "Library ID",
                    "Rank_in_Group": "Rank in Group"
                }, inplace=True)
                st.session_state.final_fragments_df = final_fragments_df

    if st.session_state.previous_runs and st.session_state.found_barcodes_df_main is not None:
        # ---------- Section: Found Barcodes Summary ----------
        with st.expander("Found Barcodes Summary", expanded=False):
            st.write("This section provides a summary of the barcodes found in the selected runs.")
            st.write("The table below shows the barcodes that were found in the library and their counts across samples.")
            found_barcodes_df = st.session_state.found_barcodes_df_main
            num_barcodes = len(found_barcodes_df)
            num_barcodes_in_lib = found_barcodes_df['Found in Library'].apply(lambda x: isinstance(x, str) and x != "Not Found").sum()

            # Create two columns
            col1, col2 = st.columns([0.2, 1])

            with col2:
                st.metric(label="Total Barcodes", value=num_barcodes)
                st.metric(label="Found in Library", value=num_barcodes_in_lib)

            with col1:
                sizes = [num_barcodes_in_lib, num_barcodes - num_barcodes_in_lib]
                labels = ['Found in Library', 'No match in Library']
                colors = ["#03254D", "#8DC8F0"]

                fig, ax = plt.subplots(figsize=(6, 6))
                wedges, texts, autotexts = ax.pie(
                    sizes,
                    labels=labels,
                    autopct='%1.1f%%',
                    startangle=90,
                    colors=colors,
                    wedgeprops=dict(width=0.15, edgecolor='w'),
                )

                for text in texts:
                    text.set_fontsize(22)
                for autotext in autotexts:
                    autotext.set_fontsize(22)

                ax.axis('equal')
                st.pyplot(fig, use_container_width=True)

            st.dataframe(found_barcodes_df, use_container_width=True)
            csv1 = found_barcodes_df.to_csv(index=False)
            st.download_button(
                label="📥 Download Found Barcodes Summary as CSV",
                data=csv1,
                file_name="found_barcodes_summary.csv",
                mime="text/csv"
            )

            # Display the Barcode Coverage per Sample
            st.subheader("Barcode Coverage per Sample")
            barcode_coverage_per_sample = st.session_state.barcode_coverage_per_sample
            st.dataframe(barcode_coverage_per_sample, use_container_width=True)
            # Download the barcode coverage per sample
            csv2 = barcode_coverage_per_sample.to_csv(index=False)
            st.download_button(
                label="📥 Download Barcode Coverage per Sample as CSV",
                data=csv2,
                file_name="barcode_coverage_per_sample.csv",
                mime="text/csv"
            )

        # ---------- Section: Find Best Fragments ----------
        with st.expander("Find Best Fragments", expanded=False):
            st.write("""
                Here we can compare fragments across different groups.  
                You can select groups that you want to include.  
                And choose if you want to select by Rank or by Reads per Million Reads.
            """)
            
            df = st.session_state.final_fragments_df

            if df is not None:
                options = df["Group"].unique()
                # sort options alphabetically descending
                options.sort()
                options = options[-1::-1] 
                col1, col2 = st.columns(2)
                groups = col1.multiselect("Groups for Comparison", options)
                metric = col2.selectbox("Select Metric", ["Rank in Group", "Reads per Million Reads"], index=0)

                if groups:
                    display_df = df[df["Group"].isin(groups)][["Peptide", "Library ID"]].drop_duplicates()

                    for group in groups:
                        gd = df[df["Group"] == group][["Peptide", "Library ID", "Reads per Million Reads", "Rank in Group"]]
                        gd = gd.rename(columns={
                            "Reads per Million Reads": f"{group}: reads per million",
                            "Rank in Group": f"Rank in {group}"
                        })
                        display_df = pd.merge(display_df, gd, on=["Peptide", "Library ID"], how="left")

                    # Keep only metric-related columns
                    if metric == "Rank in Group":
                        display_df = display_df[["Peptide", "Library ID",] + [col for col in display_df.columns if col.startswith("Rank in")]]
                        display_df = display_df.sort_values(by=[col for col in display_df.columns if col.startswith("Rank in")], ascending=True)
                    elif metric == "Reads per Million Reads":
                        display_df = display_df[["Peptide", "Library ID",] + [col for col in display_df.columns if "reads per million" in col]]
                        display_df = display_df.sort_values(by=[col for col in display_df.columns if "reads per million" in col], ascending=False)
                    # Display grid
                    gb = GridOptionsBuilder.from_dataframe(display_df)
                    gb.configure_default_column(filterable=True, sortable=True, resizable=True)
                    gridOptions = gb.build()

                    response = AgGrid(
                        display_df,
                        gridOptions=gridOptions,
                        height=400,
                        theme="alpine",
                        update_mode=GridUpdateMode.MODEL_CHANGED,
                    )

                    # Download filtered view
                    filtered_df = pd.DataFrame(response['data'])
                    csv = filtered_df.to_csv(index=False)
                    st.download_button(
                        label="📥 Download Displayed Table as CSV",
                        data=csv,
                        file_name="filtered_best_fragments.csv",
                        mime="text/csv"
                    )

    # === Project Barcode Coverage Summary ===
    st.subheader("Project Barcode Coverage Summary")
    st.write("This section gives you a quick overview about the barcode coverage of the project.")

    # Project selection
    selected_project_2 = st.selectbox(
        "Project",
        projects,
        index=projects.index(st.session_state.get('project', projects[0])),
        help="Select a project to analyze its library sequencing data.",
        key="selected_project_2"
    )

    # Dynamically set the project_dir and libraries_dir based on the selected project
    project_dir = os.path.join(PROJECTS_DIR, selected_project_2)
    libraries_dir = os.path.join(project_dir, "Libraries")

    # Load (cached) barcode summary specific to the selected project
    summary_df = load_barcode_summary(libraries_dir)

    # Display
    st.dataframe(summary_df, use_container_width=True)

    # Download
    csv = summary_df.to_csv(index=False)
    st.download_button(
        label="📥 Download Project Barcode Coverage Summary as CSV",
        data=csv,
        file_name=f"{selected_project_2}_barcode_coverage_summary.csv",
        mime="text/csv"
    )
# ------------------ Report Page ------------------

elif st.session_state.page == 'report':
    st.title(f"Library Report: `{st.session_state.project}_{st.session_state.library}_{st.session_state.run}`")
    if st.button("⬅️ Back to Selection"):
        go_back()
        # reload the page to reset state
        st.rerun()

    input_dir = st.session_state.input_dir

    # SECTION: Library Sequencing Summary
    with st.expander("Library Sequencing Summary", expanded=False):
        summary_path = get_file_path(input_dir, "sequencing_summary.csv")
        df = get_df_from_session_or_file("summary_df", summary_path)
        if df is not None:
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Sequencing Summary")
                st.table(df)
            with col2:
                pie_chart = os.path.join(input_dir, "plots/barcode_pie_chart.png")
                if os.path.exists(pie_chart):
                    st.subheader("Barcode Distribution")
                    st.image(pie_chart, use_container_width=True)
        # Display the the Logos
        barcode_reads_logo = os.path.join(input_dir, "plots/barcode_reads.svg")
        fragment_reads_logo = os.path.join(input_dir, "plots/fragment_reads.svg")
        if os.path.exists(barcode_reads_logo) and os.path.exists(fragment_reads_logo):
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Barcode Reads Logo")
                st.image(barcode_reads_logo, use_container_width=True)
            with col2:
                st.subheader("Fragment Reads Logo")
                st.image(fragment_reads_logo, use_container_width=True)
        unique_barcodes_logo = os.path.join(input_dir, "plots/unique_barcodes_logo.svg")
        unique_fragments_logo = os.path.join(input_dir, "plots/unique_fragments_logo.svg")
        if os.path.exists(unique_barcodes_logo) and os.path.exists(unique_fragments_logo):
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Unique Barcodes Logo")
                st.image(unique_barcodes_logo, use_container_width=True)
            with col2:
                st.subheader("Unique Fragments Logo")
                st.image(unique_fragments_logo, use_container_width=True)
        else:
            st.error("❌ Sequencing summary file not found.")

    # SECTION: Found Sample Barcodes Summary
    with st.expander("Summary of all Found Barcodes across Samples", expanded=False):
        report_path = get_file_path(input_dir, "found_barcode_report.csv")
        df = get_df_from_session_or_file("barcode_report_df", report_path)
        if df is not None:
            df = add_group_column(df, input_dir)
            df = df.rename(columns={
                "Name": "Sample Name",
                "BC_reads": "Found Barcode Reads",
                "unique_BC": "Found Unique Barcodes",
                "matched_BC_reads": "Barcode Reads matched to Library",
                "matched_unique_BC": "Unique Barcodes matched to Library"
            })
            gb = GridOptionsBuilder.from_dataframe(df)
            gb.configure_default_column(filterable=True, sortable=True, resizable=True)
            gb.configure_selection("single", use_checkbox=True)
            gridOptions = gb.build()
            barcode_summary = AgGrid(df, gridOptions=gridOptions, height=400, theme="alpine", fit_columns_on_grid_load=True,
                                    update_mode=GridUpdateMode.SELECTION_CHANGED)

            selected_rows = barcode_summary.get("selected_rows", [])
            if selected_rows:
                selected_sample = selected_rows[0]["Sample Name"]
                filename = get_file_path(input_dir, f"found_barcodes/unique_barcodes_{selected_sample}.csv")
                found_df = get_df_from_session_or_file(f"unique_barcodes_{selected_sample}", filename)
                if found_df is not None:
                    found_df = found_df.rename(columns={
                        "BC": "Barcode",
                        "Count": "Raw read count",
                        "RNAcount_per_mil_reads": "RNA count per mil reads",
                        "Found_in_Lib": "Found in Library Sequencing"
                    })
                    st.subheader(f"All found Barcodes in {selected_sample}")
                    st.dataframe(found_df, use_container_width=True)
        else:
            st.warning("❌ `found_barcode_report.csv` not found.")

    # ---------- Section: All Found Barcodes ----------
    with st.expander("List of most abundant Barcodes", expanded=False):
        st.write("This section provides a summary of the most abundant barcodes found across all samples. It includes all barcodes with a mean read count per million reads greater than 1.\n")
        found_barcodes_path = get_file_path(input_dir, "found_barcodes.csv")
        report_path = get_file_path(input_dir, "found_barcode_report.csv")

        found_df = get_df_from_session_or_file("found_barcodes_df", found_barcodes_path)
        report_df = get_df_from_session_or_file("barcode_report_df", report_path)

        if "annotated_report_df" not in st.session_state and report_df is not None:
            st.session_state["annotated_report_df"] = add_group_column(report_df, input_dir)

        annotation = st.session_state.get("annotated_report_df", None)

        if found_df is not None and annotation is not None:
            num_samples = len(found_df.columns) - 5
            found_df = found_df[found_df['Count_per_mil_reads_mean'] > 1]

            rename_map = {
                "BC": "Found Barcode",
                "Found_in_Lib": "Found in Library",
                "Count_per_mil_reads_mean": "Mean Reads per Million Reads",
                "Max_per_mil_reads": "Max Reads per Million Reads in a Sample",
                "Samples_Found_In": f"Samples Found in ... out of {num_samples}"
            }

            for col_name in found_df.columns:
                if col_name in annotation['Name'].values:
                    matched_group = annotation.loc[annotation['Name'] == col_name, 'Group']
                    if not matched_group.empty:
                        new_name = f"{matched_group.values[0]} ({col_name})"
                        rename_map[col_name] = new_name

            found_df.rename(columns=rename_map, inplace=True)

            num_barcodes = len(found_df)
            num_barcodes_in_lib = found_df['Found in Library'].sum()
            # Create three columns
            col1, col2 = st.columns([0.2, 1])

            # Left metric
            with col2:
                st.metric(label="Total Barcodes", value=num_barcodes)
                st.metric(label="Found in Library", value=num_barcodes_in_lib)

            # Center donut chart
            with col1:
                sizes = [num_barcodes_in_lib, num_barcodes - num_barcodes_in_lib]
                labels = ['Found in Library', 'No match in Library']
                colors = ["#03254D", "#8DC8F0"]  # Dark blue and light blue

                fig, ax = plt.subplots(figsize=(5, 5))  # Smaller chart
                wedges, texts, autotexts = ax.pie(
                    sizes,
                    labels=labels,
                    autopct='%1.1f%%',
                    startangle=90,
                    colors=colors,
                    wedgeprops=dict(width=0.25, edgecolor='w'),
                )
                # Set font size for labels
                for text in texts:
                    text.set_fontsize(22)

                # Set font size for autopct values
                for autotext in autotexts:
                    autotext.set_fontsize(22)
                ax.axis('equal')
                st.pyplot(fig, use_container_width=True) # Adjusted size for better fit

            gb = GridOptionsBuilder.from_dataframe(found_df)
            gb.configure_default_column(filterable=True, sortable=True, resizable=True)
            AgGrid(found_df, gridOptions=gb.build(), height=400, theme="alpine", allow_unsafe_jscode=True)
        else:
            st.warning("❌ Required barcode files not found.")

    # ---------- Section: Final Fragment Summary ----------
    with st.expander("Summary of all Fragments that could be traced with Barcodes", expanded=False):
        st.write("""
        This section provides a summary of all fragments that were detected by matching barcodes found in the samples to the library barcodes.\n
        You can take a closer look at the fragments by selecting them from the table below.\n
        If *1* Fragment is selected, you will see the details of the fragment across all groups.\n
        If *2* Fragments are selected, you will see a comparison of the two fragments across all groups.\n
        """)
        summary_path = get_file_path(input_dir, "final_fragments_summary.csv")
        df = get_df_from_session_or_file("final_fragments_df", summary_path)

        if df is not None:
            # Filter only once
            if "filtered_fragments_df" not in st.session_state:
                df_filtered = df[df["Group"] != "Plasmid_Library"]
                df_filtered = df_filtered[[
                    "Group", "Peptide", "LLinker", "Sequence", "RLinker",
                    "RNAcount", "RNAcount_ratio", "BC_count", "BC", "in_subsets"
                ]]
                df_filtered.rename(columns={
                    "RNAcount": "Reads per Million Reads",
                    "RNAcount_ratio": "Abundance in Group",
                    "BC_count": "Barcode Count",
                    "BC": "Barcode(s)",
                    "in_subsets": "Fraction of groups represented in subsets"
                }, inplace=True)
                st.session_state["filtered_fragments_df"] = df_filtered

            df_filtered = st.session_state["filtered_fragments_df"]

            gb = GridOptionsBuilder.from_dataframe(df_filtered)
            gb.configure_default_column(filterable=True, sortable=True, resizable=True)
            gb.configure_column("Group", header_name="Group", filter="agSetColumnFilter")
            gb.configure_selection("multiple", use_checkbox=True)
            gridOptions = gb.build()

            final_fragments = AgGrid(
                df_filtered,
                gridOptions=gridOptions,
                height=600,
                theme="alpine",
                update_mode=GridUpdateMode.SELECTION_CHANGED,
                allow_unsafe_jscode=True,
                fit_columns_on_grid_load=True,
            )

            selected_rows = final_fragments.get("selected_rows", [])

            if selected_rows:
                selected_df = pd.DataFrame(selected_rows)
                if len(selected_df) > 2:
                    st.warning("⚠️ Only the first 2 rows will be used.")
                    selected_df = selected_df.head(2)

                if len(selected_df) == 1:
                    peptide = selected_df["Peptide"].iloc[0]
                    match = df_filtered[df_filtered["Peptide"] == peptide]
                    st.subheader(f"Details for Peptide: {peptide}")
                    st.dataframe(match, use_container_width=True)

                elif len(selected_df) == 2:
                    p1, p2 = selected_df["Peptide"].iloc[0], selected_df["Peptide"].iloc[1]
                    m1 = df_filtered[df_filtered["Peptide"] == p1]
                    m2 = df_filtered[df_filtered["Peptide"] == p2]
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
                    st.dataframe(
                        cmp[["Group", f"Reads per Million Reads ({p1})", f"Reads per Million Reads ({p2})", "Ratio"]],
                        use_container_width=True
                    )
        else:
            st.warning("❌ `final_fragments_summary.csv` not found.")

            
    # ---------- Section: Find Best Fragments ----------
    with st.expander("Find Best Fragments", expanded=False):
        st.write("""
            Here we can compare fragments across different groups.  
            You can select groups that you want to include.  
            And choose if you want to select by Rank or by Reads per Million Reads.
        """)

        # Use cached dataframe
        df = st.session_state.get("final_fragments_df", None)

        if df is not None:
            # Rename columns once
            if "best_fragments_df" not in st.session_state:
                df = df.rename(columns={
                    "in_subsets": "Fraction of groups represented in subsets",
                    "RNAcount": "Mean Reads per Million Reads",
                    "Rank_in_Group": "Rank in Group"
                })
                st.session_state["best_fragments_df"] = df
            else:
                df = st.session_state["best_fragments_df"]

            options = df["Group"].unique()
            # sort options alphabetically descending
            options.sort()
            options = options[-1::-1] 
            col1, col2 = st.columns(2)
            groups = col1.multiselect("Groups for Comparison", options)
            metric = col2.selectbox("Select Metric", ["Rank in Group", "Mean Reads per Million Reads"])

            if groups:
                display_df = df[df["Group"].isin(groups)][["Peptide"]].drop_duplicates()

                for group in groups:
                    gd = df[df["Group"] == group][["Peptide", "Mean Reads per Million Reads", "Rank in Group"]]
                    gd = gd.rename(columns={
                        "Mean Reads per Million Reads": f"{group}: reads per million",
                        "Rank in Group": f"Rank in {group}"
                    })
                    display_df = pd.merge(display_df, gd, on="Peptide", how="left")

                # Keep only metric-related columns
                if metric == "Rank in Group":
                    display_df = display_df[["Peptide"] + [col for col in display_df.columns if col.startswith("Rank in")]]
                elif metric == "Mean Reads per Million Reads":
                    display_df = display_df[["Peptide"] + [col for col in display_df.columns if "reads per million" in col]]

                # Display grid
                gb = GridOptionsBuilder.from_dataframe(display_df)
                gb.configure_default_column(filterable=True, sortable=True, resizable=True)
                gridOptions = gb.build()

                response = AgGrid(
                    display_df,
                    gridOptions=gridOptions,
                    height=400,
                    theme="alpine",
                    update_mode=GridUpdateMode.MODEL_CHANGED,
                )

                # Download filtered view
                filtered_df = pd.DataFrame(response['data'])
                csv = filtered_df.to_csv(index=False)
                st.download_button(
                    label="📥 Download Displayed Table as CSV",
                    data=csv,
                    file_name="filtered_best_fragments.csv",
                    mime="text/csv"
                )
        else:
            st.warning("❌ `final_fragments_summary.csv` not found.")

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
            st.warning("❌ Heatmaps not found.")