#!/bin/bash

###############################################################################
# Script Name: combine_undetermined_with_libraries.sh
# Description: This script combines the undetermined reads with the reads from
#              each library in the specified folder. The script requires the
#              path to the folder containing the libraries and the undetermined
#              reads as input.
#
# Usage:
#   ./combine_undetermined_with_libraries.sh /path/to/folder
#
# Parameters:
#   - /path/to/folder: (Required) Path to the folder containing the libraries
#
# Author: Jaro Steindorff
###############################################################################

# Check if the folder path is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 /path/to/folder"
    exit 1
fi

# Define the folder path
folder_path="$1"

# Check if the required directories exist
if [ ! -d "$folder_path" ]; then
    echo "Error: The folder '$folder_path' does not exist."
    exit 1
fi

if [ ! -d "$folder_path/undetermined" ]; then
    echo "Error: The 'undetermined' folder is missing in '$folder_path'."
    exit 1
fi

# Define libraries and file types
libraries=("p005" "p006" "p007")
file_types=("barcode" "fragment")

# Loop through each library and file type
for library in "${libraries[@]}"; do
    for file_type in "${file_types[@]}"; do
        # Construct input and output file paths
        input1="${folder_path}${library}/barcode_fragment/${file_type}_${library}.fastq.gz"
        input2="${folder_path}undetermined/barcode_fragment/${file_type}_undetermined_${library}.fastq.gz"
        output="${folder_path}${library}/barcode_fragment/combined_${file_type}_${library}.fastq.gz"

        # Check if input files exist
        if [ ! -f "$input1" ]; then
            echo "Warning: Input file '$input1' not found. Skipping."
            continue
        fi

        if [ ! -f "$input2" ]; then
            echo "Warning: Input file '$input2' not found. Skipping."
            continue
        fi

        # Concatenate the files
        cat "$input1" "$input2" > "$output"
        
        # Notify user of the concatenation
        echo "Processed ${output}"
    done
done

echo "All files processed successfully!"
