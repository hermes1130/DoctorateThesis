library("DEXSeq")
library("stringr")
library("stringi")
library("GenomicFeatures")
library("GenomicRanges")

pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
list.files(pythonScriptsDir)
cat(pythonScriptsDir) #/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/DEXSeq/python_scripts

#Make a list of input files
inDir <- "data/DEX_input/Neurons"
countFiles <- list.files(inDir, pattern=".txt$", full.names=TRUE)
countFiles

#Define the DEX gtf file
flattenedFile <- list.files("data/", pattern="_DEX.gtf$", full.names=TRUE)
flattenedFile

#prepare two sample tables

x <- c(0,2,3,5,10,15,25)
Samples_h <- c()
Samples_g <- c()
for (i in 1:3){
  g_rep <- stri_paste("Gorilla", i, "_day", x)
  Samples_g <- c(Samples_g, g_rep)

  h_rep <- stri_paste("Human", i, "_day", x)
  Samples_h <- c(Samples_h, h_rep)
  
  Samples <- c(Samples_g, Samples_h)
}

sampleTable <- data.frame(
row.names = Samples,
species = c( rep("Gorilla", 21),
rep("Human", 21)))

sampleTable_t <- data.frame(
row.names = Samples,
species = c( rep("Gorilla", 21),
rep("Human", 21)), 
timepoint = c(rep(x, 3)))

#construct an DEXSeqDataSet object 
dxd <- DEXSeqDataSetFromHTSeq(
    countFiles,
    sampleData=sampleTable,
    design= ~ sample + exon + species:exon,
    flattenedfile=flattenedFile )
    

colData(dxd)

dxd[which(stringr::str_detect(rownames(dxd), "ENSG00000169554")),] #110 exon bins for ZEB2

head( counts(dxd), 5 )
#Notice that the number of columns is 14, the first seven (we have seven samples) corresponding to the
#number of reads mapping to out exonic regions and the last seven correspond to the sum of the counts
#mapping to the rest of the exons from the same gene on each sample.

head( featureCounts(dxd), 5 )
#We can also access only the first five rows from the count belonging to the exonic regions 

head( rowRanges(dxd), 5 )
#this table contains information on the annotation data, such as gene and exon IDs, genomic
#coordinates of the exon, and the list of transcripts that contain an exon.

sampleAnnotation( dxd )
#the sample annotation

dxd <- estimateSizeFactors( dxd )
#DEXSeq uses the same method as DESeq and DESeq2

dxd <- estimateDispersions( dxd )
#Dispersion estimation - takes around 20-30 min
#png("output/Figures/2024_11_13_Neurons_DispEstsPlot.png", width = 1000, height = 1500)
#svg("output/Figures/2024_11_13_Neurons_DispEstsPlot.svg", width = 10, height = 15)
plotDispEsts( dxd )
#dev.off()
#per-exon dispersion estimates versus the mean normalised count, 
#the resulting fitted valuesand the a posteriori (shrinked) dispersion estimates 

dxd <- testForDEU( dxd )
#Testing for differential exon usage

dxd_Neuron <- estimateExonFoldChanges( dxd, fitExpToVar="species")

saveRDS(dxd_Neuron, file = "output/preliminary/dxd_Neuron.rds")
dxd_Neuron <- readRDS("output/preliminary/dxd_Neuron.rds")

# Ensure the species column is a factor and reorder it as desired
dxd_Neuron$species <- factor(dxd_Neuron$species, levels = c("Human", "Gorilla")) 

rst3 <- DEXSeqResults( dxd_Neuron )
rst3

table ( rst3$padj < 0.05 )
#how many exonic regions are significant with a false discovery rate of 5%
# FALSE   TRUE 
#134100  16993 
table ( tapply( rst3$padj < 0.05, rst3$groupID, any ))
#how many genes are affected
# FALSE   TRUE 
# 1791  4911 

#png("output/Figures/2024_11_14_Neurons_DEU.png", width = 1000, height = 1000)
#svg("output/Figures/2024_11_14_Neurons_DEU.svg", width = 10, height = 15)
plotDEXSeq(rst3, "ENSG00000169554.23", fitExpToVar="species", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
#dev.off()
#png("output/Figures/2024_11_14_Neurons_DEU_tranascripts.png", width = 1000, height = 1000)
#svg("output/Figures/2024_11_14_Neurons_DEU_tranascripts.svg", width = 10, height = 15)
plotDEXSeq(rst3, "ENSG00000169554.23", fitExpToVar="species", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
#This gene contains more than 40 transcripts annotated, only the first 40 will be plotted
#dev.off()

###Which bins belong to which exons?
# Extract information on bins from the DEXSeqResults object
bins_info_3 <- as.data.frame(rst3)
gene_bins_3 <- bins_info_3[bins_info_3$groupID == "ENSG00000169554.23", ]
nrow(gene_bins_3) #110
table (gene_bins_3$padj < 0.05 ) #nothing is significant
#FALSE  TRUE 
#   45     1 

# Assuming `gene_bins_3` is a data frame, create a GRanges object for it
bins_gr_3 <- GRanges(
  seqnames = gene_bins_3$genomicData.seqnames,
  ranges = IRanges(start = gene_bins_3$genomicData.start, end = gene_bins_3$genomicData.end),
  strand = gene_bins_3$genomicData.strand
)
seqlevels(bins_gr_3)

# Load GTF file
txdb <- makeTxDbFromGFF("data/gencode.v46.chr_patch_hapl_scaff.annotation.gtf", format = "gtf")
# Extract exons for ZEB2
#zeb2_exons <- exons(txdb, filter = list(gene_id = "ENSG00000169554.23")) 

# Merge overlapping exons to get unique ranges
#zeb2_exons_unique <- reduce(zeb2_exons)

# Filter by the main transcript ID, if known (replace with actual transcript ID)
main_transcript_exons <- exons(txdb, filter = list(tx_name = "ENST00000627532.3"))

# Check the resulting exon count
main_transcript_exons
seqlevels(main_transcript_exons)
# Find overlaps between bins and ZEB2 exons
overlaps_t_3 <- findOverlaps(bins_gr_3, main_transcript_exons, type = "any")

# Add exon information to `gene_bins_2` using overlaps
gene_bins_3$exonID <- rep(NA, nrow(gene_bins_3))  # Placeholder for exonID

# Annotate the bins with exon IDs where overlaps are found
gene_bins_3$exonID[queryHits(overlaps_t_3)] <- subjectHits(overlaps_t_3)

unique(gene_bins_3$exonID)

sum(is.na(gene_bins_3$exonID)) #67 bins out of 110 bins are intergenic regions

# Convert 'transcripts' column to character by collapsing each list into a string
gene_bins_3$transcripts <- sapply(gene_bins_3$transcripts, function(x) paste(x, collapse = ","))
write.table(gene_bins_3, file = "output/DEXresults_Neurons.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

