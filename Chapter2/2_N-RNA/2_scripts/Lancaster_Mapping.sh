#!/bin/bash
#SBATCH -J Lancaster2021_mapping
#SBATCH -D /data/scratch/hermes1130
#SBATCH -o Lancaster2021_mapping.%j.out
#SBATCH --partition=big
#SBATCH --cpus-per-task=5
#SBATCH --nodes=2
#SBATCH --mem=5000
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=hermes1130@zedat.fu-berlin.de

#conda activation
source activate RNA_seq_2024

#INPUT DATA PREP

# put the unique sample names into an array
sample_files_h=($(cut -f 1 /home/hermes1130/Lancaster/sample_name_Lancaster2021.txt))

# Replace the base path with the correct /cutadapt path
for i in "${!sample_files_h[@]}"; do
    sample_files_h[$i]="${sample_files_h[$i]/\/data\/scratch\/hermes1130\/input\/Lancaster2021\//\/data\/scratch\/hermes1130\/input\/Lancaster2021\/cutadapt\/}"
done

# print out all the element in the array
echo "Human samples: ${sample_files_h[@]}"

# Align human samples with Hisat2 for single-end RNA-seq
for file in "${sample_files_h[@]}"
do
    R1="${file}_Atrimmed.fastq.gz"  # Single-end file
    echo "Aligning human sample: $file"
    hisat2 -q --score-min L,0,-0.6 --rdg 7,5 -x /data/scratch/hermes1130/ref/human/GRCh38_OnlyKnwonChr -U "$R1" -S "${file}_human.sam"
done

