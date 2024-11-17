library("DEXSeq")
library("stringr")
library("GenomicFeatures")
library("GenomicRanges")

pythonScriptsDir = system.file( "python_scripts", package="DEXSeq" )
list.files(pythonScriptsDir)
cat(pythonScriptsDir) #/Library/Frameworks/R.framework/Versions/4.2-arm64/Resources/library/DEXSeq/python_scripts

#Make a list of input files
inDir <- "data/DEX_input/NCC"
countFiles <- list.files(inDir, pattern=".txt$", full.names=TRUE)
countFiles

#Define the DEX gtf file
flattenedFile <- list.files("data/", pattern="_DEX.gtf$", full.names=TRUE)
flattenedFile

#prepare a sample table
sampleTable <- data.frame(
row.names = c( "Human1_rep1", "Human1_rep2", 
"Human2_rep1", "Human2_rep2", 
"Human3_rep1", "Human3_rep2", 
"Chimp1_rep1", "Chimp1_rep2",
"Chimp2_rep1", "Chimp2_rep2" ),
species = c("Human", "Human", "Human","Human", "Human", "Human", 
"Chimp",  "Chimp", "Chimp", "Chimp"))

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
#png("output/Figures/2024_11_10_NCC_DispEstsPlot.png", width = 1000, height = 1500)
#svg("output/Figures/2024_11_10_NCC_DispEstsPlot.svg", width = 10, height = 15)
plotDispEsts( dxd )
#dev.off()
#per-exon dispersion estimates versus the mean normalised count, 
#the resulting fitted valuesand the a posteriori (shrinked) dispersion estimates 

dxd <- testForDEU( dxd )
#Testing for differential exon usage

dxd_NCC <- estimateExonFoldChanges( dxd, fitExpToVar="species")

saveRDS(dxd_NCC, file = "output/preliminary/dxd_NCC.rds")
dxd_NCC <- readRDS("output/preliminary/dxd_NCC.rds")

# Ensure the species column is a factor and reorder it as desired
dxd_NCC$species <- factor(dxd_NCC$species, levels = c("Human", "Chimp")) 

rst2 <- DEXSeqResults( dxd_NCC )
rst2
rst2[which(stringr::str_detect(rownames(rst2), "ENSG00000169554")),] 

table ( rst2$padj < 0.05 )
#how many exonic regions are significant with a false discovery rate of 5%
# FALSE   TRUE 
#130198    406 
table ( tapply( rst2$padj < 0.05, rst2$groupID, any ))
#how many genes are affected
#FALSE  TRUE 
#397   278 

#png("output/Figures/2024_11_10_NCC_DEU.png", width = 1000, height = 1000)
#svg("output/Figures/2024_11_10_NCC_DEU.svg", width = 10, height = 15)
plotDEXSeq(rst2, "ENSG00000169554.23", fitExpToVar="species", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
#dev.off()
#png("output/Figures/2024_11_10_NCC_DEU_tranascripts.png", width = 1000, height = 1000)
#svg("output/Figures/2024_11_10_NCC_DEU_tranascripts.svg", width = 10, height = 15)
plotDEXSeq(rst2, "ENSG00000169554.23", fitExpToVar="species", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
#This gene contains more than 40 transcripts annotated, only the first 40 will be plotted
#dev.off()

###Which bins belong to which exons?
# Extract information on bins from the DEXSeqResults object
bins_info_2 <- as.data.frame(rst2)
gene_bins_2 <- bins_info_2[bins_info_2$groupID == "ENSG00000169554.23", ]
nrow(gene_bins_2)
table (gene_bins_2$padj < 0.05 ) #nothing is significant
#FALSE 
#    2 

# Assuming `gene_bins_2` is a data frame, create a GRanges object for it
bins_gr_2 <- GRanges(
  seqnames = gene_bins_2$genomicData.seqnames,
  ranges = IRanges(start = gene_bins_2$genomicData.start, end = gene_bins_2$genomicData.end),
  strand = gene_bins_2$genomicData.strand
)
seqlevels(bins_gr_2)

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
overlaps_t_2 <- findOverlaps(bins_gr_2, main_transcript_exons, type = "any")

# Add exon information to `gene_bins_2` using overlaps
gene_bins_2$exonID <- rep(NA, nrow(gene_bins_2))  # Placeholder for exonID

# Annotate the bins with exon IDs where overlaps are found
gene_bins_2$exonID[queryHits(overlaps_t_2)] <- subjectHits(overlaps_t_2)

unique(gene_bins_2$exonID)

sum(is.na(gene_bins_2$exonID)) #67 bins out of 110 bins are intergenic regions

# Convert 'transcripts' column to character by collapsing each list into a string
gene_bins_2$transcripts <- sapply(gene_bins_2$transcripts, function(x) paste(x, collapse = ","))
write.table(gene_bins_2, file = "output/DEXresults_NCC.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
