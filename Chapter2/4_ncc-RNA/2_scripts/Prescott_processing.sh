#!/bin/bash
#SBATCH -J Prescott2015_processing
#SBATCH -D /data/scratch/hermes1130
#SBATCH -o Prescott2015_processing.%j.out
#SBATCH --partition=big
#SBATCH --cpus-per-task=5
#SBATCH --nodes=2
#SBATCH --mem=5000
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=hermes1130@zedat.fu-berlin.de

#conda activation
source activate RNA_seq_2024

#convert sam to bam
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*.sam | parallel --progress --eta -j 9 'samtools view -S -b {} > {}.bam'
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*.sam.bam | rename 's/.sam.bam/.bam/'

#Mapped + MAPQ > 10 (can go higher if you want)
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*.bam | parallel --progress --eta -j 9 'samtools view -bq 10 {} > {.}_mapQ10.bam'

# Sort Bam Files
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*_mapQ10.bam | parallel --progress --eta -j 9 'samtools sort {} > {.}_sorted.bam'
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*_mapQ10_sorted.bam | rename 's/_mapQ10_sorted.bam/_sorted.bam/'

#Namesorting
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*_sorted.bam | parallel --progress --eta -j 9 'samtools sort -n {} > {.}_namesorted.bam'
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*_sorted_namesorted.bam | rename 's/_sorted_namesorted.bam/_namesorted.bam/'

#Fixmating
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*_namesorted.bam | parallel --progress --eta -j 9 'samtools fixmate -rc {} {.}_fixed.bam'
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*_namesorted_fixed.bam | rename 's/_namesorted_fixed.bam/_fixed.bam/'

#Re-sorting
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*_fixed.bam | parallel --progress --eta -j 9 'samtools sort {} >  {.}_resorted.bam'
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*_fixed_resorted.bam | rename 's/_fixed_resorted.bam/_resorted.bam/'

#Indexing
ls /data/scratch/hermes1130/input/Prescott2015/mapping/*_resorted.bam | parallel --progress --eta -j 9 'samtools index {}'
