""" This file contains the  config dictionary that is used to store the configuration parameters for the pipeline. """

data_dir = "0_data"

config_S1 = {
    "input_file": "input/DNA-lib_RetrogradeTransport.fasta",
    "wSet": "input/wSet.csv",
    "14_aa_overhangs": ["AACCTCCAGAGAGGCAACGCT", "GCCAGACAAGCAGCTACCGCA"],
    "14_aa_G4S_overhangs": ["AACCTCCAGAGAGGCAACGGAGGCGGAGGAAGT", "GGAGGCGGCGGAAGCAGACAAGCAGCTACCGCA"],
    "14_aa_A5_overhangs": ["AACCTCCAGAGAGGCAACGCTGCTGCAGCAGCC", "GCAGCTGCAGCTGCCAGACAAGCAGCTACCGCA"],
    "22_aa_overhangs": ["AACCTCCAGAGAGGCAACGCT", "GCCAGACAAGCAGCTACCGCA"],
    "output_name": data_dir + "/SortedFragments.txt",
}

config_S2 = {
    "in_name_P5": data_dir + "/fastq_files/DNA_pscAAVlib_1.fastq.gz",
    "in_name_P7": data_dir + "/fastq_files/DNA_pscAAVlib_2.fastq.gz",
    "out_dir": data_dir + "/barcode_fragment",
    "out_name": "DNA_pscAAVlib_1",
    "bbduk2_args_BC" : [        
        "overwrite=true",
        "k=18",
        "mink=18",
        "hammingdistance=2",
        "findbestmatch=t",
        "rcomp=f",
        "findbestmatch=f",
        "qhdist=1",
        "minavgquality=0",
        "maxns=0",
        "minlength=18",
        "maxlength=22",
        "lliteral=GGCCTAGCGGCCGCTTTACTT",
        "rliteral=ATAACTTCGTATAATGTATGC",
    ],
    "bbduk2_args_Frag" : [
        "overwrite=true",
        "k=18",
        "mink=18",
        "rcomp=f",
        "qhdist=1",
        "maskmiddle=t",
        "hammingdistance=2",
        "findbestmatch=t",
        "minlength=38",
        "maxlength=78",
        "ordered=t",
        "lliteral=AGCAACCTCCAGAGAGGCAACG",
        "rliteral=CAGACAAGCAGCTACCGCAGAT",
    ],
}

config_S3 = {
    "original_seq_file": config_S1["input_file"],
    "wSet": "input/wSet.csv",
    "fragment_seq_file": config_S1["output_name"],
    "constitutive_backbone_sequences": ["aacctccagagaggcaacg", "cagacaagcagctaccgca"],
    "linker_dict": {
        "14aaG4S": "GAGGCGGAGGAAGT",
        "14aaA5": "CTGCTGCAGCAGCC",
    },
    "length_dict": {
        "14aa": 46,
        "22aa": 70,
    },
    "trim_dict": {
        "14aa": [2,44],
        "14aaG4S": [14,56],
        "14aaA5": [14,56],
        "22aa": [2,68],
    },
    "out_name_LUT": data_dir + "/LUT.csv",
    "out_name": data_dir + "/fragment_pos.csv",
}

config_S4 = {
    "in_name_LUT": config_S3["out_name_LUT"],
    "barcode_file": config_S2["out_dir"] + "/barcode_" + config_S2["out_name"] + ".fastq.gz",
    "fragment_file": config_S2["out_dir"] + "/fragment_" + config_S2["out_name"] + ".fastq.gz",
    "out_name": data_dir + "/library_fragments.csv",
}

config_S5 = {
    "input_table": config_S4["out_name"],
    "fragments_pos": config_S3["out_name"],
    "in_name_LUT": config_S3["out_name_LUT"],
    "sample_inputs": "input/load_list.csv",
    "sample_directory": data_dir + "/fastq_files",
    "log_file_path": data_dir + "/found_barcode_report.csv",
    "output_dir": data_dir + "/found_barcodes",
    "bbduk2_args" : [        
        "overwrite=true",
        "k=12",
        "mink=12",
        "hammingdistance=2",
        "findbestmatch=t",
        "trd=t",
        "rcomp=f",
        "skipr2=t",
        "findbestmatch=f",
        "qhdist=0",
        "minavgquality=0",
        "ordered=t",
        "maxns=0",
        "minlength=18",
        "maxlength=22",
        "lliteral=GGCCTAGCGGCCGCTTTACTT",
        "rliteral=ATAACTTCGTATA"
    ],
}

config_S6 = {
    "input_table": config_S4["out_name"],
    "fragments_pos": config_S3["out_name"],
    "output_table": data_dir + "/annotated_library_fragments.csv",
}

config_S7 = {
    "original_seq_file": config_S1["input_file"],
    "input_dir": config_S5["output_dir"],
    "library_fragments": config_S6["output_table"],
    "subsets": {
        "Infective_AAVs": ['exclude','DNA_AAVlib_DNAse_30cpc', 'DNA_AAVlib_DNAse_3cpc','DNA_pscAAVlib', 'DNA_pscAAVlib_Prep2'],
        "DNAse_resistant_AAVs": ['include', 'DNA_AAVlib_DNAse_30cpc','DNA_AAVlib_DNAse_3cpc'],
    },
    "trim_dict": config_S3["trim_dict"],
    "output_table": data_dir + "/final_fragments_summary.csv",
}

# create a lookup dictionary for the configuration dictionaries
config_lookup = {
    "S1": config_S1,
    "S2": config_S2,
    "S3": config_S3,
    "S4": config_S4,
    "S5": config_S5,
    "S6": config_S6,
    "S7": config_S7,
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