""" This file contains the config dictionary that are used to store the configuration parameters for the pipeline. """

# Define the data directory where the input and output files are stored
data_dir = "raav-60/undetermined"
# Define the directory where the logs are stored
log_dir = data_dir + "/logs/"
# Define the length of the barcode and fragment sequences in DNA bases
bc_len = 27
frag_len = 21

config_S2 = {
    # input file names for the P5 and P7 fastq files P5 is the barcode and P7 is the fragment
    "in_name_barcode": data_dir + "/fastq_files/Undetermined_R1.fastq.gz",
    "in_name_fragment": data_dir + "/fastq_files/Undetermined_R2.fastq.gz",
    "input_file": data_dir + "/input/reference_seq_p005.fasta",
    # output directory and name for the barcode and fragment files once they have been extracted
    "out_dir": data_dir + "/barcode_fragment",
    "out_name": "undetermined_p005",
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
        "lliteral=GTACGTCTGAACTTGGGACT",
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

# create a lookup dictionary for the configuration dictionaries
config_lookup = {
    "S2": config_S2,
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
