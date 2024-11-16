#!/bin/bash

# Create the species name
for file in *_AA.fasta; do
  species_name="${file%_AA.fasta}"
  echo "$species_name"
  # Use sed to replace the first line with the species name (adding '>' for FASTA format)
  sed -i '' "1s/.*/>${species_name}/" "$file"
done

# Create an output file
output_file="combined_AA.fasta"

# Append the contents of each FASTA file to the output file
for file in *.fasta; do
  cat "$file" >> "$output_file"
done

# Define the order in an array
order=(
  "Homo_sapiens"
  "Pan_troglodytes" "Pan_paniscus"
  "Gorilla_gorilla_gorilla"
  "Pongo_abelii" "Pongo_pygmaeus"
  "Nomascus_leucogenys" "Symphalangus_syndactylus"
  "Hylobates_moloch"
  "Macaca_mulatta" "Macaca_fascicularis"
  "Macaca_thibetana_thibetana" "Macaca_nemestrina"
  "Papio_anubis" "Theropithecus_gelada"
  "Cercocebus_atys"
  "Chlorocebus_sabaeus"
  "Rhinopithecus_bieti" "Rhinopithecus_roxellana"
  "Trachypithecus_francoisi"
  "Piliocolobus_tephrosceles"
  "Colobus_angolensis_palliatus"
  "Cebus_imitator" "Sapajus_apella"
  "Saimiri_boliviensis_boliviensis"
  "Callithrix_jacchus" "Aotus_nancymaae"
  "Carlito_syrichta" "Nycticebus_coucang"
  "Propithecus_coquereli" "Microcebus_murinus" "Lemur_catta"
  "Otolemur_garnettii" "Mus_musculus"
)

# Read through the FASTA file and store sequences by header
declare -A sequences
current_header=""

while read -r line; do
  if [[ $line == ">"* ]]; then
    # Store the header without '>'
    current_header="${line#>}"
    sequences["$current_header"]="$line"$'\n'
  else
    # Append sequence lines to the current header's sequence
    sequences["$current_header"]+="$line"$'\n'
  fi
done < combined_AA.fasta

# Output in the desired order
for species in "${order[@]}"; do
  if [[ -n ${sequences[$species]} ]]; then
    echo "${sequences[$species]}"
  else
    echo "Warning: $species not found in file."
  fi
done > ordered_combined_AA.fasta
