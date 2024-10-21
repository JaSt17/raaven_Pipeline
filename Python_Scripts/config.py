""" This file contains the  config dictionary that is used to store the configuration parameters for the pipeline. """

config_S1 = {
    "input_file": "input/DNA-lib_RetrogradeTransport.fasta",
    "wSet": "input/wSet.csv",
    "14_aa_overhangs": ["AACCTCCAGAGAGGCAACGCT", "GCCAGACAAGCAGCTACCGCA"],
    "14_aa_G4S_overhangs": ["AACCTCCAGAGAGGCAACGGAGGCGGAGGAAGT", "GGAGGCGGCGGAAGCAGACAAGCAGCTACCGCA"],
    "14_aa_A5_overhangs": ["AACCTCCAGAGAGGCAACGCTGCTGCAGCAGCC", "GCAGCTGCAGCTGCCAGACAAGCAGCTACCGCA"],
    "22_aa_overhangs": ["AACCTCCAGAGAGGCAACGCT", "GCCAGACAAGCAGCTACCGCA"],
    "output_name": "0_data/SortedFragments_all_python.txt",
}

config_S2 = {
    "in_name_P5": "0_data/fastq_files/DNA_pscAAVlib_1.fastq.gz",
    "in_name_P7": "0_data/fastq_files/DNA_pscAAVlib_2.fastq.gz",
    "out_dir": "0_data/barcode_fragment",
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
    "fragment_seq_file": "0_data/SortedFragments_all_python.txt",
    "constitutive_backbone_sequences": ["aacctccagagaggcaacg", "cagacaagcagctaccgca"],
    "linker_dict": {
        "14aaG4S": "GAGGCGGAGGAAGT",
        "14aaA5": "CTGCTGCAGCAGCC",
    },
    "length_dict": {
        "14aa": 14,
        "22aa": 22,
    },
    "trim_dict": {
        "14aa": [2,44],
        "14aaG4S": [14,56],
        "14aaA5": [14,56],
        "22aa": [2,68],
    },
    "out_name": "0_data/alignedLibraries.csv",
}

# create a lookup dictionary for the configuration dictionaries
config_lookup = {
    "S1": config_S1,
    "S2": config_S2,
    "S3": config_S3,
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