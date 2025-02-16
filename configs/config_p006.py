""" This file contains the config dictionary that are used to store the configuration parameters for the pipeline. """

# Define the data directory where the input and output files are stored
data_dir = "raav-60/p006"
# Define the directory where the logs are stored
log_dir = data_dir + "/logs/"
# Define the length of the barcode and fragment sequences in DNA bases
bc_len = 27
frag_len = 21

config_S2 = {
    # input file names for the P5 and P7 fastq files P5 is the barcode and P7 is the fragment
    "in_name_barcode": data_dir + "/fastq_files/p006_R1.fastq.gz",
    "in_name_fragment": data_dir + "/fastq_files/p006_R2.fastq.gz",
    # output directory and name for the barcode and fragment files once they have been extracted
    "out_dir": data_dir + "/barcode_fragment",
    "out_name": "p006",
    # arguments for the bbduk2 tool to extract the barcode and fragment sequences
    "bbduk2_args_BC" : [
        "k=20",
        "hammingdistance=1",
        "overwrite=true",
        "findbestmatch=t",
        "rcomp=f",
        "minavgquality=0",
        "maxns=0",
        f"minlength={bc_len}",
        f"maxlength={bc_len}",
        "ordered=t",
        "lliteral=GCCATCCTCTTATCTCGTGG",
        "rliteral=ATAACTTCGTATAATGTATG",
    ],
    "bbduk2_args_Frag" : [
        "k=20",
        "hammingdistance=1",
        "overwrite=true",
        "findbestmatch=t",
        "maskmiddle=t",
        "rcomp=f",
        "minavgquality=0",
        "maxns=0",
        f"minlength={frag_len}",
        f"maxlength={frag_len}",
        "ordered=t",
        "lliteral=ACCTCCAGAGAGGCAACGCT",
        "rliteral=GCCAGACAAGCAGCTACCGC",
    ],
    "log_dir": log_dir,
}

config_S3 = {
    # input file names are extracted from the previous step
    "barcode_file": config_S2["out_dir"] + "/barcode_" + config_S2["out_name"] + ".fastq.gz",
    "fragment_file": config_S2["out_dir"] + "/fragment_" + config_S2["out_name"] + ".fastq.gz",
    # Do we want to allwo single read barcodes
    "single_read": True,
    # Do we want to allow chimeric barcodes
    "chimeric_read": False,
    # threshold for the ratio of the most frequent barcode to all found barcodes for chimeric barcode detection
    "threshold": 0.9,
    # the chunk size determains how many sequences are read in at once and can be set to a smaller number if memory is an issue
    "chunk_size": 20000000,
    # output file name for the library barcodes
    "out_name": data_dir + "/library_barcodes.csv",
    "log_dir": log_dir,
}

config_S4 = {
    # input file names are extracted from the previous step
    "input_table": config_S3["out_name"],
    "in_name_LUT": None,
    "chunk_size": config_S3["chunk_size"],
    "bc_len": bc_len,
    "db": data_dir + "/barcode_db.fasta",
    # input csv file containing the file names of all samples that should be used for barcode extraction
    "sample_inputs": data_dir + "/input/load_list.csv",
    # directory containing the fastq files for the samples
    "sample_directory": "raav-60/sample_fastq",
    # filename for the log file that will be created and show how many barcodes were found in each sample
    "log_file_path": data_dir + "/found_barcode_report.csv",
    # output directory for the found barcodes csv files
    "output_dir": data_dir + "/found_barcodes",
    # arguments for the bbduk2 tool to extract the barcodes from the samples
    "bbduk2_args" : [        
        "k=20",
        "hammingdistance=1",
        "overwrite=true",
        "findbestmatch=t",
        "maskmiddle=t",
        "rcomp=f",
        "minavgquality=0",
        "maxns=0",
        f"minlength={bc_len}",
        f"maxlength={bc_len}",
        "ordered=t",
        "lliteral=GCCATCCTCTTATCTCGTGG",
        "rliteral=ATAACTTCGTATAATGTATG",
    ],
    "log_dir": log_dir,
}

config_S5 = {
    # input file names are extracted from the previous step
    "input_table": config_S3["out_name"],
    "in_name_LUT": None,
    # output file name for the library barcodes with their information form the LUT
    "output_table": data_dir + "/pos_library_barcodes.csv",
    "log_dir": log_dir,
}

config_S6 = {
    # input file names are extracted from the previous step
    "input_dir": config_S4["output_dir"],
    "sample_inputs": config_S4["sample_inputs"],
    "library_fragments": config_S5["output_table"],
    # group name for the library
    "library_name": "Plasmid_Library",
    # dictionary containing the information about the different subsets that should be created
    # the key is the name of the subset and the value is a list of the fragments that should be included
    "subsets": {
        "Infective_AAVs": ['exclude','DNAse_resistant_AAVs','Plasmid_Library'],
    },
    # output file name for the final fragments summary
    "output_table": data_dir + "/final_fragments_summary.csv",
    "log_dir": log_dir,
}

# create a lookup dictionary for the configuration dictionaries
config_lookup = {
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