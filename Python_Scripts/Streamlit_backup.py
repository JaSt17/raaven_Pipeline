# ---------- Section: Fragment Group Comparison ----------
with st.expander("Compare Fragment Abundance Across Samples/Tissues", expanded=False):
    st.write("""
        Here we can compare the effectiveness of fragments across different groups.\n
        You can select a group to compare against other groups. The first group will be shown in the first column, while the other groups will be shown in the following columns.\n
        The comparison is based on the reads per million reads for each fragment in the selected groups.\n
        """)
    summary_file = os.path.join(input_dir, "final_fragments_summary.csv")
    if os.path.exists(summary_file):
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
    else:
        st.warning("❌ `final_fragments_summary.csv` not found.")