#!/bin/bash
#SBATCH -J Prescott2015_mapping
#SBATCH -D /data/scratch/hermes1130
#SBATCH -o Prescott2015_mapping.%j.out
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

#define the metadata input
meta_file="/data/scratch/hermes1130/input/Prescott2015/metadata.txt"

#define output files for human and chimp sample names
human_sample_file="/data/scratch/hermes1130/Prescott_sample_name_human.txt"
chimp_sample_file="/data/scratch/hermes1130/Prescott_sample_name_chimp.txt"

#define the directory path to prepend
input_dir="/data/scratch/hermes1130/input/Prescott2015/cutadapt/"

#clear the output files if they already exist
> "$human_sample_file"
> "$chimp_sample_file"

# Debugging: Ensure file variables are correctly set
echo "Human samples file: $human_sample_file"
echo "Chimp samples file: $chimp_sample_file"

#read the metadata file line by line
while IFS=$'\t' read -r study_accession sample_accession experiment_accession run_accession tax_id scientific_name fastq_ftp submitted_ftp sra_ftp bam_ftp; do
  # Check the species and append the sample name to the appropriate file
	if [ "$scientific_name" == "Homo sapiens" ]; then
		echo "Writing ${input_dir}${run_accession} to $human_sample_file"
		echo "${input_dir}${run_accession}" >> "$human_sample_file"
	elif [ "$scientific_name" == "Pan troglodytes" ]; then
		echo "Writing ${input_dir}${run_accession} to $chimp_sample_file"
		echo "${input_dir}${run_accession}" >> "$chimp_sample_file"
	fi
done < "$meta_file"

echo "Sample lists created: $human_sample_file and $chimp_sample_file"

# put the unique sample names into an array
sample_files_h=($(cut -f 1 /data/scratch/hermes1130/Prescott_sample_name_human.txt))
sample_files_c=($(cut -f 1 /data/scratch/hermes1130/Prescott_sample_name_chimp.txt))

# print out all the element in the array
echo "Human samples: ${sample_files_h[@]}"
echo "Chimp samples: ${sample_files_c[@]}"

# Align human samples with Hisat2
for file in "${sample_files_h[@]}"
do
    R1="${file}_Atrimmed_1.fastq.gz"
    R2="${file}_Atrimmed_2.fastq.gz"
    echo "Aligning human sample: $file"
    hisat2 -q --score-min L,0,-0.6 --rdg 7,5 -x /data/scratch/hermes1130/ref/human/GRCh38_OnlyKnwonChr -1 "$R1" -2 "$R2" -S "${file}_human.sam"
done

# Align chimp samples with Hisat2
for file in "${sample_files_c[@]}"
do
    R1="${file}_Atrimmed_1.fastq.gz"
    R2="${file}_Atrimmed_2.fastq.gz"
    echo "Aligning chimp sample: $file"
    hisat2 -q --score-min L,0,-0.6 --rdg 7,5 -x /data/scratch/hermes1130/ref/chimp/panTro6.dna -1 "$R1" -2 "$R2" -S "${file}_chimp.sam"
done

#convert sam to bam
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*.sam | parallel --progress --eta -j 9 'samtools view -S -b {} > {}.bam'
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*.sam.bam | rename 's/.sam.bam/.bam/'

#Mapped + MAPQ > 10 (can go higher if you want)
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*.bam | parallel --progress --eta -j 9 'samtools view -bq 10 {} > {.}_mapQ10.bam'

# Sort Bam Files
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*_mapQ10.bam | parallel --progress --eta -j 9 'samtools sort {} > {.}_sorted.bam'
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*_mapQ10_sorted.bam | rename 's/_mapQ10_sorted.bam/_sorted.bam/'

#Namesorting
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*_sorted.bam | parallel --progress --eta -j 9 'samtools sort -n {} > {.}_namesorted.bam'
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*_sorted_namesorted.bam | rename 's/_sorted_namesorted.bam/_namesorted.bam/'

#Fixmating
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*_namesorted.bam | parallel --progress --eta -j 9 'samtools fixmate -rc {} {.}_fixed.bam'
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*_namesorted_fixed.bam | rename 's/_namesorted_fixed.bam/_fixed.bam/'

#Re-sorting
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*_fixed.bam | parallel --progress --eta -j 9 'samtools sort {} >  {.}_resorted.bam'
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*_fixed_resorted.bam | rename 's/_fixed_resorted.bam/_resorted.bam/'

#Indexing
#ls /data/scratch/hermes1130/rawdata_RNA/chimp/*_resorted.bam | parallel --progress --eta -j 9 'samtools index {}'


