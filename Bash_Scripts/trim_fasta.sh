#!/bin/bash

# Usage: ./trim_fasta.sh input.fasta > output.fasta

input_file=$1

awk '
BEGIN {seq = ""} 
/^>/ { 
    if (seq != "") {
        print substr(seq, 3, length(seq) - 4);
    }
    print $0; 
    seq = "";
    next;
} 
{ seq = seq $0 } 
END { 
    if (seq != "") {
        print substr(seq, 3, length(seq) - 4);
    }
}' "$input_file"
