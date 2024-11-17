#!/bin/bash
#SBATCH -J Lancaster2021_processing
#SBATCH -D /data/scratch/hermes1130
#SBATCH -o Lancaster2021_processingV.%j.out
#SBATCH --partition=big
#SBATCH --cpus-per-task=5
#SBATCH --nodes=2
#SBATCH --mem=5000
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=hermes1130@zedat.fu-berlin.de

#conda activation
source activate RNA_seq_2024

# Define the directory
input_dir="/data/scratch/hermes1130/input/Lancaster2021/mapping/"

# Convert SAM to BAM
find "$input_dir" -name "*.sam" | parallel --progress --eta -j 9 'samtools view -S -b {} > {.}.bam'

# Filter BAM files by MAPQ > 10
find "$input_dir" -name "*.bam" | parallel --progress --eta -j 9 'samtools view -bq 10 {} > {.}_mapQ10.bam'

# Sort BAM files
find "$input_dir" -name "*_mapQ10.bam" | parallel --progress --eta -j 9 'samtools sort {} -o {.}_sorted.bam'

# Name sorting of BAM files
find "$input_dir" -name "*_sorted.bam" | parallel --progress --eta -j 9 'samtools sort -n {} -o {.}_namesorted.bam'

# Fixmate
find "$input_dir" -name "*_namesorted.bam" | parallel --progress --eta -j 9 'samtools fixmate -rc {} {.}_fixed.bam'

# Re-sorting BAM files
find "$input_dir" -name "*_fixed.bam" | parallel --progress --eta -j 9 'samtools sort {} -o {.}_resorted.bam'

# Indexing BAM files
find "$input_dir" -name "*_resorted.bam" | parallel --progress --eta -j 9 'samtools index {}'

echo "Pipeline completed."