#!/bin/bash

# Define libraries and file types
libraries=("p005" "p006" "p007")
file_types=("barcode" "fragment")

# Loop through each library and file type
for library in "${libraries[@]}"; do
    for file_type in "${file_types[@]}"; do
        # Construct input and output file paths
        input1="${library}/${file_type}_fragment/small_${file_type}_${library}.fastq.gz"
        input2="undetermined/${file_type}_fragment/${file_type}_undetermined_${library}.fastq.gz"
        output="${library}/${file_type}_fragment/${file_type}_${library}.fastq.gz"

        # Concatenate the files
        cat "$input1" "$input2" > "$output"
        
        # Notify user of the concatenation
        echo "Processed ${output}"
    done
done

echo "All files processed successfully!"
