#!/bin/bash

# Check if input file is provided
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <input_bed_file>"
  exit 1
fi

# Load the tabix module
module load tabixpp

input_file=$1
sorted_file="${input_file}.sorted"
gz_file="${sorted_file}.gz"

# Sort the BED file
sort -k1,1 -k2,2n "$input_file" > "$sorted_file"

# Compress the sorted file
bgzip -f "$sorted_file"

# Index the compressed file
tabix -f -p bed "$gz_file"

echo "Processing complete: $gz_file and index created."


