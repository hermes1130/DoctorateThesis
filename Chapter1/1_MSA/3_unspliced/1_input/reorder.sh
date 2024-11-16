#!/bin/bash

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

# Directory for temporary files
temp_dir=$(mktemp -d)

# Split the FASTA file by headers into separate files in temp_dir
awk '/^>/ {filename=sprintf("%s/%s.fasta", "'$temp_dir'", substr($0, 2)); print > filename; next} {print >> filename}' combined.fasta

# Create the ordered output file
output_file="ordered_combined.fasta"
> "$output_file"  # Empty the output file if it exists

# Concatenate each species' file in the specified order
for species in "${order[@]}"; do
  species_file="$temp_dir/$species.fasta"
  if [[ -f "$species_file" ]]; then
    cat "$species_file" >> "$output_file"
  else
    echo "Warning: $species not found in file."
  fi
done

# Clean up temporary files
rm -rf "$temp_dir"

echo "Reordering complete. Output saved to $output_file."
