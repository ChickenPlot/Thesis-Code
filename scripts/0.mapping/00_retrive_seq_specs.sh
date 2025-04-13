#!/bin/bash
fq_dir="../X204SC22090037-Z01-F002/01.RawData"
output_file="../X204SC22090037-Z01-F002/01.RawData/sequencers_specs.txt"

# Create or overwrite the output file
> "$output_file"

# Loop through each .fq file in the directory
for fq_file in "$fq_dir"/*.fq; do
    # Extract the first line and append it to the output file
    filename=$(basename "$fq_file")
    echo "$filename" >> "$output_file"
    head -n 1 "$fq_file" >> "$output_file"
done