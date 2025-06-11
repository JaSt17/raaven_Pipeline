#------------------------------------------------------------------------------------
# Define the data & save directory where the input and output files are stored
data_dir = "Projects/Bluejay/Seq_Data"
save_dir = "Projects/Bluejay/Libraries"
# Name of the library and run
library_name = "p034"
run = "Plasmid_run_1"
# Define the length of the barcode and fragment sequences in DNA bases
bc_len = 27
frag_len = 27
linker_length = 3  # Lenght of the Linker left and right to the fragment sequence

# Define the number of possible fragments that can be created from the library
num_possible_frag = 180000  # This is an arbitrary number, adjust as needed

# Define the literals for the barcode and fragment sequences
barcode_left_literal = "CTATCGAGTG" # Unique !!!
barcode_right_literal = "ATAACTTCGT"
fragment_left_literal = "GAGTGCCCAA"
fragment_right_literal = "GCACAGGCGC"

# Settings for Libary read usage:
# Should single read barcodes be used?
single_read_barcodes = True
# Should chimeric barcodes be used?
chimeric_barcodes = False
# Should starcode reduction be used?
starcode_reduction = False
# threshold for the retrival of chimeric barcodes
threshold = 1
#------------------------------------------------------------------------------------

# configuration for Step 1 in the pipeline
config_S1 = {
    # input file containing the DNA sequences to create the library from
    "input_file": save_dir + f"/{library_name}/{run}/input/reference_seq.fasta",
    # dictionary containing the information about the different structures with
    # their name as the key and then the length, frequency
    "structure_dict": {
        "7aa": {"length": 7, "freq": 1}},
    # Library ID for the library so we can combine multiple libraries in the future
    "LibID": library_name,
    # output file names for the LUT csv and the list of all inserted fragments
    "output_csv": save_dir + f"/{library_name}/{run}/intermediate_files/LUT.csv",
    "output_name": save_dir + f"/{library_name}/{run}/intermediate_files/SortedFragments.txt",
    "log_dir": save_dir + f"/{library_name}/{run}/logs/",
}

config_S2 = {
    # input file names for the P5 and P7 fastq files P5 is the barcode and P7 is the fragment
    "in_name_barcode": data_dir + f"/{run}/{library_name}_R1.fastq.gz",
    "in_name_fragment": data_dir + f"/{run}/{library_name}_R2.fastq.gz",
    "input_file": config_S1["input_file"],
    "draw_sequence_logos": True,  # Whether to draw sequence logos for the barcodes and fragments
    # output directory and name for the barcode and fragment files once they have been extracted
    "out_dir": save_dir + f"/{library_name}/{run}/barcode_fragment",
    "out_name": library_name,
    # arguments for the bbduk2 tool to extract the barcode and fragment sequences
    "bbduk2_args_BC" : [
        "k=10",
        "hammingdistance=1",
        "overwrite=true",
        "findbestmatch=t",
        "rcomp=f",
        "minavgquality=0",
        "maxns=0",
        f"minlength={bc_len}",
        f"maxlength={bc_len}",
        "ordered=t",
        f"lliteral={barcode_left_literal}",
        f"rliteral={barcode_right_literal}",
    ],
    "bbduk2_args_Frag" : [
        "k=10",
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
        f"lliteral={fragment_left_literal}",
        f"rliteral={fragment_right_literal}",
    ],
    "log_dir": save_dir + f"/{library_name}/{run}/logs/",
}

config_S3 = {
    # input file names are extracted from the previous step
    "in_name_LUT": config_S1["output_csv"],
    "barcode_file": config_S2["out_dir"] + "/barcode_" + config_S2["out_name"] + ".fastq.gz",
    "fragment_file": config_S2["out_dir"] + "/fragment_" + config_S2["out_name"] + ".fastq.gz",
    "library_name": library_name,
    # Do we want to allwo single read barcodes
    "single_read": single_read_barcodes,
    # Do we want to allow chimeric barcodes
    "chimeric_read": chimeric_barcodes,
    # Do we want to use starcode reduction
    "starcode": starcode_reduction,
    # threshold for the ratio of the most frequent barcode to all found barcodes for chimeric barcode detection
    "threshold": threshold,
    # the chunk size determains how many sequences are read in at once and can be set to a smaller number if memory is an issue
    "chunk_size": 20000000,
    # output file name for the library barcodes
    "out_name": save_dir + f"/{library_name}/{run}/intermediate_files/library_barcodes.csv",
    "log_dir": save_dir + f"/{library_name}/{run}/logs/",
}

config_S4 = {
    # input file names are extracted from the previous step
    "input_table": config_S3["out_name"],
    "in_name_LUT": None,
    "chunk_size": config_S3["chunk_size"],
    "bc_len": bc_len,
    "starcode": config_S3["starcode"],
    "db": save_dir + f"/{library_name}/{run}/intermediate_files/barcode_db.fasta",
    # input csv file containing the file names of all samples that should be used for barcode extraction
    "sample_inputs": data_dir + "/Samples/annotation.csv",
    # directory containing the fastq files for the samples
    "sample_directory": data_dir + "/Samples",
    # filename for the log file that will be created and show how many barcodes were found in each sample
    "log_file_path": save_dir + f"/{library_name}/{run}/found_barcode_report.csv",
    # output directory for the found barcodes csv files
    "output_dir": save_dir + f"/{library_name}/{run}/found_barcodes",
    # arguments for the bbduk2 tool to extract the barcodes from the samples
    "bbduk2_args" : [
        "k=10",
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
        f"lliteral={barcode_left_literal}",
        f"rliteral={barcode_right_literal}",
    ],
    "log_dir": save_dir + f"/{library_name}/{run}/logs/",
}

config_S5 = {
    # input file names are extracted from the previous step
    "input_table": config_S3["out_name"],
    "in_name_LUT": config_S1["output_csv"],
    # output file name for the library barcodes with their information form the LUT
    "output_table": save_dir + f"/{library_name}/{run}/intermediate_files/pos_library_barcodes.csv",
    "log_dir": save_dir + f"/{library_name}/{run}/logs/",
}

config_S6 = {
    # input file names are extracted from the previous step
    "original_seq_file": config_S1["input_file"],
    "input_dir": config_S4["output_dir"],
    "sample_inputs": config_S4["sample_inputs"],
    "library_fragments": config_S5["output_table"],
    "linker_length": linker_length,
    "plot_dir": save_dir + f"/{library_name}/{run}/plots",
    "array_size": num_possible_frag,  # size of the array for the summary plots
    # group name for the library
    "library_name": "Plasmid_Library",
    # dictionary containing the information about the different subsets that should be created
    # the key is the name of the subset and the value is a list of the fragments that should be included
    "subsets": {
        "Infective_AAVs": ['exclude','DNAse_resistant_AAVs','Plasmid_Library'],
    },
    # output file name for the final fragments summary
    "output_table": save_dir + f"/{library_name}/{run}/final_fragments_summary.csv",
    "log_dir": save_dir + f"/{library_name}/{run}/logs/",
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