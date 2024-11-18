#!/bin/bash
#SBATCH -J Prescott2015_QC
#SBATCH -D /data/scratch/hermes1130
#SBATCH -o Prescott2015_QC.%j.out
#SBATCH --partition=big
#SBATCH --cpus-per-task=5
#SBATCH --nodes=2
#SBATCH --mem=50000
#SBATCH --time=24:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=hermes1130@zedat.fu-berlin.de

#conda activation
source activate RNA_seq_2024

#run the initial run with fastqc
#fastqc -o /data/scratch/hermes1130/input/Prescott2015/fastqc -f fastq /data/scratch/hermes1130/input/Prescott2015/*.fastq.gz

#make a txt file containing all the sample names
#ls /data/scratch/hermes1130/input/Prescott2015/*.fastq.gz | sed "s/_1.fastq.gz//;s/_2.fastq.gz//;" | awk '!seen[$0]++' > sample_name_Prescott2015.txt

# put the unique sample names into an array
#sample_files=($(cut -f 1 /home/hermes1130/sample_name_Prescott2015.txt))
# print out all the element in the array
#echo "${sample_files[@]}"

#Quality trimming is done BEFORE any adapter trimming
#remove the reads with lower quality of 20 and smaller read length of 21
#for file in "${sample_files[@]}"
#do
#	R1="${file}"_1.fastq.gz
#	R2="${file}"_2.fastq.gz
#	cutadapt -q 20 -m 21 --pair-filter=any -o ${file}_QC_1.fastq.gz -p ${file}_QC_2.fastq.gz "$R1" "$R2"
#done


#remove the adapter sequences collected from the initial fastqc run
#for file in "${sample_files[@]}"
#do
#        R1="${file}"_QC_1.fastq.gz
#        R2="${file}"_QC_2.fastq.gz
#        cutadapt -a file:/data/scratch/hermes1130/input/Prescott2015/fastqc/Prescott_adapter_sequences.fasta -A file:/data/scratch/hermes1130/input/Prescott2015/fastqc/Prescott_adapter_sequences.fasta -o ${file}_Atrimmed_1.fastq.gz -p ${file}_Atrimmed_2.fastq.gz "$R1" "$R2"
#done

#move the output files to the corresponding folder
#mkdir /data/scratch/hermes1130/input/Prescott2015/cutadapt
#mv *QC* /data/scratch/hermes1130/input/Prescott2015/cutadapt
#mv *Atrimmed* /data/scratch/hermes1130/input/Prescott2015/cutadapt

#run the fastqc for the checkup
fastqc -o /data/scratch/hermes1130/input/Prescott2015/fastqc -f fastq /data/scratch/hermes1130/input/Prescott2015/cutadapt/*Atrimmed*.fastq.gz
