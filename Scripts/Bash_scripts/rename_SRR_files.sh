#!/bin/bash

# Define the path to the CSV file
csv_file="$1"

# Define the directory where FASTQ files are located
fastq_dir="$2"

# Read the CSV file line by line
while IFS=, read -r srr library_name; do
    # Skip the header line
    if [[ "$srr" == "Run" ]]; then
        continue
    fi

    # Ensure the FASTQ directory exists
    if [ ! -d "$fastq_dir" ]; then
        echo "FASTQ directory not found: $fastq_dir"
        exit 1
    fi

    # Loop through FASTQ files that match the SRR number pattern
    for suffix in _1 _2; do
        file="${fastq_dir}/${srr}${suffix}.fastq"
        # Check if the file exists
        if [ -f "$file" ]; then
            # Extract the file extension
            extension="${file##*.}"
            # Construct the new file name
            new_name="${fastq_dir}/${library_name}${suffix}.${extension}"
            # Rename the file
            mv "$file" "$new_name"
            echo "Renamed $file to $new_name"
        else
            echo "File not found for SRR $srr: $file"
        fi
    done

done < "$csv_file"