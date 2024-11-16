#!/bin/bash

# Create the species name
for file in *.fasta; do
  species_name="${file%.fasta}"
  echo "$species_name"
  # Use sed to replace the first line with the species name (adding '>' for FASTA format)
  sed -i '' "1s/.*/>${species_name}/" "$file"
done

# Create an output file
output_file="combined.fasta"

# Append the contents of each FASTA file to the output file
for file in *.fasta; do
  cat "$file" >> "$output_file"
done