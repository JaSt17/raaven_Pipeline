#!/bin/bash
set -euo pipefail

###############################################################################
# Script Name: combine_undetermined_with_libraries.sh
# Description: Combines undetermined reads with library reads for analysis.
# Author: Jaro Steindorff
###############################################################################

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/folder"
    exit 1
fi

folder_path="$1"

if [ ! -d "$folder_path" ]; then
    echo "Error: The folder '$folder_path' does not exist."
    exit 1
fi

if [ ! -d "$folder_path/Undetermined" ]; then
    echo "Error: The 'Undetermined' folder is missing in '$folder_path'."
    exit 1
fi

# set the libraries and plasmid run
libraries=("p005" "p006" "p007")
plasmid_run="Plasmid_run_1"
file_types=("barcode" "fragment")

for library in "${libraries[@]}"; do
    for file_type in "${file_types[@]}"; do
        input1="${folder_path}/${library}/${plasmid_run}/barcode_fragment/${file_type}_${library}.fastq.gz"
        input2="${folder_path}/Undetermined/${plasmid_run}/barcode_fragment/${file_type}_Undetermined_${library}.fastq.gz"
        output="${folder_path}/${library}/${plasmid_run}/barcode_fragment/combined_${file_type}_${library}.fastq.gz"

        if [ ! -f "$input1" ]; then
            echo "Warning: Input file '$input1' not found. Skipping."
            continue
        fi

        if [ ! -f "$input2" ]; then
            echo "Warning: Input file '$input2' not found. Skipping."
            continue
        fi

        cat "$input1" "$input2" > "$output"

        mv "$input1" "${folder_path}/${library}/${plasmid_run}/barcode_fragment/small_${file_type}_${library}.fastq.gz"
        mv "$output" "$input1"

        echo "Processed ${input1} and renamed original to small_${file_type}_${library}.fastq.gz"
    done
done

echo "✅ All files processed successfully!"
