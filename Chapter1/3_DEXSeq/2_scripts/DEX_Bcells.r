library("DEXSeq")
library("stringr")
library("GenomicFeatures")
library("GenomicRanges")

pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
list.files(pythonScriptsDir)
cat(pythonScriptsDir) #/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/DEXSeq/python_scripts

#Make a list of input files
inDir <- "data/DEX_input/B-cells"
countFiles <- list.files(inDir, pattern="_ctrl.txt$", full.names=TRUE)
countFiles

#Define the DEX gtf file
flattenedFile <- list.files("data/", pattern="_DEX.gtf$", full.names=TRUE)
flattenedFile

#prepare a sample table
sampleTable <- data.frame(
row.names = c( "Human_GM18558_ctrl", "Human_GM18960_ctrl", "Human_GM19240_ctrl",
"Chimp_Judith_ctrl", "Chimp_Leo_ctrl", "Chimp_Maryke_ctrl", 
"Orang_Guchi_ctrl", "Orang_Jaqo_ctrl", "Orang_JingJing_ctrl" ),
species = c("Human", "Human", "Human",
"Chimp", "Chimp", "Chimp", 
"Orang",  "Orang", "Orang"))

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
#png("output/Figures/2024_11_08_Bcell_DispEstsPlot.png", width = 1000, height = 1500)
#svg("output/Figures/2024_11_08_Bcell_DispEstsPlot.svg", width = 10, height = 15)
plotDispEsts( dxd )
#dev.off()
#per-exon dispersion estimates versus the mean normalised count, 
#the resulting fitted valuesand the a posteriori (shrinked) dispersion estimates 

dxd <- testForDEU( dxd )
#Testing for differential exon usage

dxd_B <- estimateExonFoldChanges( dxd, fitExpToVar="species")

saveRDS(dxd_B, file = "output/preliminary/dxd_B.rds")
dxd_B <- readRDS("output/preliminary/dxd_B.rds")

# Ensure the species column is a factor and reorder it as desired
dxd_B$species <- factor(dxd_B$species, levels = c("Human", "Chimp", "Orang")) 

rst1 <- DEXSeqResults( dxd_B )
rst1
rst1[which(stringr::str_detect(rownames(rst1), "ENSG00000169554")),] 

table ( rst1$padj < 0.05 )
#how many exonic regions are significant with a false discovery rate of 5%
# FALSE   TRUE 
#261901  28610 
table ( tapply( rst1$padj < 0.05, rst1$groupID, any ))
#how many genes are affected
#FALSE  TRUE 
# 1101  6526                    

#png("output/Figures/2024_11_08_Bcell_DEU.png", width = 1000, height = 1000)
#svg("output/Figures/2024_11_08_Bcell_DEU.svg", width = 10, height = 15)
plotDEXSeq(rst1, "ENSG00000169554.23", fitExpToVar="species", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
#dev.off()
#png("output/Figures/2024_11_10_Bcell_DEU_tranascripts.png", width = 1000, height = 1000)
#svg("output/Figures/2024_11_10_Bcell_DEU_tranascripts.svg", width = 10, height = 15)
plotDEXSeq(rst1, "ENSG00000169554.23", fitExpToVar="species", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
#This gene contains more than 40 transcripts annotated, only the first 40 will be plotted
#dev.off()

###Which bins belong to which exons?
# Extract information on bins from the DEXSeqResults object
bins_info <- as.data.frame(rst1)
gene_bins <- bins_info[bins_info$groupID == "ENSG00000169554.23", ]
nrow(gene_bins)
table (gene_bins$padj < 0.05 ) #nothing is significant

# Assuming `gene_bins` is a data frame, create a GRanges object for it
bins_gr <- GRanges(
  seqnames = gene_bins$genomicData.seqnames,
  ranges = IRanges(start = gene_bins$genomicData.start, end = gene_bins$genomicData.end),
  strand = gene_bins$genomicData.strand
)
seqlevels(bins_gr)

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
overlaps <- findOverlaps(bins_gr, main_transcript_exons)
overlaps_t <- findOverlaps(bins_gr, main_transcript_exons, type = "any")

# Add exon information to `gene_bins` using overlaps
gene_bins$exonID <- rep(NA, nrow(gene_bins))  # Placeholder for exonID

# Annotate the bins with exon IDs where overlaps are found
gene_bins$exonID[queryHits(overlaps_t)] <- subjectHits(overlaps_t)

unique(gene_bins$exonID)

sum(is.na(gene_bins$exonID)) #67 bins out of 110 bins are intergenic regions

# Convert 'transcripts' column to character by collapsing each list into a string
gene_bins$transcripts <- sapply(gene_bins$transcripts, function(x) paste(x, collapse = ","))
write.table(gene_bins, file = "output/DEXresults_Bcells.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)