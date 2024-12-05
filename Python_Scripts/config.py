""" This file contains the config dictionary that are used to store the configuration parameters for the pipeline. """

# Define the data directory where the input and output files are stored
data_dir = "raav-60/p005"
# Define the directory where the logs are stored
log_dir = data_dir + "/logs/"

# configuration for Step 1 in the pipeline
config_S1 = {
    # input file containing the DNA sequences to create the library from
    "input_file": data_dir + "/input/reference_seq.fasta",
    # wSet file containing the hsa codon usage table
    "wSet": data_dir + "/input/wSet.csv",
    # dictionary containing the information about the different structures with
    # their name as the key and then the length, frequency, and overhangs as the values
    "structure_dict": {
        "7aa": {"length": 7, "freq": 1,
                "overhangs": ["AACCTCCAGAGAGGCAACGCT", "GCCAGACAAGCAGCTACCGCA"]}},
    # output file names for the LUT csv and the list of all inserted fragments
    "output_csv": data_dir + "/LUT.csv",
    "output_name": data_dir + "/SortedFragments.txt",
    "log_dir": log_dir,
}

config_S2 = {
    # input file names for the P5 and P7 fastq files P5 is the barcode and P7 is the fragment
    "in_name_barcode": data_dir + "/fastq_files/p005_R1.fastq.gz",
    "in_name_fragment": data_dir + "/fastq_files/p005_R2.fastq.gz",
    # output directory and name for the barcode and fragment files once they have been extracted
    "out_dir": data_dir + "/barcode_fragment",
    "out_name": "p005",
    # arguments for the bbduk2 tool to extract the barcode and fragment sequences
    "bbduk2_args_BC" : [
        "k=20",
        "hammingdistance=2",
        "overwrite=true",
        "findbestmatch=t",
        "rcomp=f",
        "qhdist=1",
        "minavgquality=0",
        "maxns=0",
        "minlength=27",
        "maxlength=27",
        "ordered=t",
        "lliteral=GTACGTCTGAACTTGGGACT",
        "rliteral=ATAACTTCGTATAATGTATGC",
    ],
    "bbduk2_args_Frag" : [
        "k=18",
        "hammingdistance=2",
        "overwrite=true",
        "findbestmatch=t",
        "maskmiddle=t",
        "rcomp=f",
        "qhdist=1",
        "minavgquality=0",
        "maxns=0",
        "minlength=25",
        "maxlength=25",
        "ordered=t",
        "lliteral=ACCAACCTCCAGAGAGGCAACG",
        "rliteral=CAGACAAGCAGCTACCGCAGAT",
    ],
    "log_dir": log_dir,
}

config_S3 = {
    # input file names are extracted from the previous step
    "in_name_LUT": config_S1["output_csv"],
    "barcode_file": config_S2["out_dir"] + "/barcode_" + config_S2["out_name"] + ".fastq.gz",
    "fragment_file": config_S2["out_dir"] + "/fragment_" + config_S2["out_name"] + ".fastq.gz",
    # threshold for the ratio of the most frequent barcode to all found barcodes for chimeric barcode detection
    "threshold": 1.0,
    # the chunk size determains how many sequences are read in at once and can be set to a smaller number if memory is an issue
    "chunk_size": 10000000,
    # output file name for the library barcodes
    "out_name": data_dir + "/library_barcodes.csv",
    "log_dir": log_dir,
}

config_S4 = {
    # input file names are extracted from the previous step
    "input_table": config_S3["out_name"],
    "in_name_LUT": config_S1["output_csv"],
    # input csv file containing the file names of all samples that should be used for barcode extraction
    "sample_inputs": data_dir + "/input/load_list.csv",
    # directory containing the fastq files for the samples
    "sample_directory": data_dir + "/fastq_files",
    # filename for the log file that will be created and show how many barcodes were found in each sample
    "log_file_path": data_dir + "/found_barcode_report.csv",
    # output directory for the found barcodes csv files
    "output_dir": data_dir + "/found_barcodes",
    # arguments for the bbduk2 tool to extract the barcodes from the samples
    "bbduk2_args" : [        
        "k=13",
        "hammingdistance=2",
        "overwrite=true",
        "findbestmatch=t",
        "rcomp=f",
        "qhdist=1",
        "trd=t",
        "skipr2=t",
        "minavgquality=0",
        "maxns=0",
        "minlength=18",
        "maxlength=22",
        "ordered=t",
        "lliteral=GGCCTAGCGGCCGCTTTACTT",
        "rliteral=ATAACTTCGTATA",
    ],
    "log_dir": log_dir,
}

config_S5 = {
    # input file names are extracted from the previous step
    "input_table": config_S3["out_name"],
    "in_name_LUT": config_S1["output_csv"],
    # output file name for the library barcodes with their information form the LUT
    "output_table": data_dir + "/pos_library_barcodes.csv",
    "log_dir": log_dir,
}

config_S6 = {
    # input file names are extracted from the previous step
    "original_seq_file": config_S1["input_file"],
    "input_dir": config_S4["output_dir"],
    "library_fragments": config_S5["output_table"],
    # group name for the library
    "library_name": "Plasmid_Library",
    # dictionary containing the information about the different subsets that should be created
    # the key is the name of the subset and the value is a list of the fragments that should be included
    "subsets": {
        "Infective_AAVs": ['exclude','DNA_AAVlib_DNAse_30cpc_1', 'DNA_AAVlib_DNAse_3cpc_1','Plasmid_Library', 'DNA_pscAAVlib_Prep2_1'],
        "DNAse_resistant_AAVs": ['include', 'DNA_AAVlib_DNAse_30cpc_1','DNA_AAVlib_DNAse_3cpc_1'],
        "Transported_AAVs": ['contains_include', "mRNA_30cpc_SN", "mRNA_30cpc_Th", "mRNA_30cpc_Ctx", "mRNA_3cpc_SN", "mRNA_3cpc_Th", "mRNA_3cpc_Ctx"],
    },
    "backbone_seq": ["aacctccagagaggcaacg", "cagacaagcagctaccgca"],
    # dictionary containing the information about how different structures should be trimmed
    "trim_dict": {
        "14aa": [2,44],
        "14aaG4S": [14,56],
        "14aaA5": [14,56],
        "22aa": [2,68],
    },
    # output file name for the final fragments summary
    "output_table": data_dir + "/final_fragments_summary.csv",
    "log_dir": log_dir,
}

# create a lookup dictionary for the configuration dictionaries
config_lookup = {
    "S1": config_S1,
    "S2": config_S2,
    "S3": config_S3,
    "S4": config_S4,
    "S5": config_S5,
    "S6": config_S6,
}

def get_config(step: str) -> dict:
    """
    Returns the configuration dictionary for the specified step.

    :param step: The step for which to return the configuration dictionary
    :return: The configuration dictionary for the specified step
    """
    try:
        return config_lookup[step]
    except KeyError:
        raise ValueError(f"Step {step} not found in the configuration lookup dictionary")