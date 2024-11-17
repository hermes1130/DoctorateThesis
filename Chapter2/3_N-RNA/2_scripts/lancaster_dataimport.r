library(here)
library(dplyr)
library(stringi)

library(DESeq2)
library(tibble)
library(DEGreport)

library(stringr)
library(ggplot2)
library(tidyr)
library(plyr)
library(webr)
#the aim of this project: getting a list of correlated genes with ZEB2 expression in human and gorilla.


#############################################################################
#Data prep
#############################################################################
#read counts preparation
human_MOF <- read.table(here("data/readcounts/GENCODE_incUnknown_genecounts_hg38_MOF_Lancaster"), header = T)
dim(human_MOF) #67062 27
gorilla_MOF <- read.table(here("data/readcounts/GENCODE_incUnknown_genecounts_KamilahGGOv0_crsmapped_MOf_Lancaster"), header = T)
dim(gorilla_MOF) #67062 27

###version number drop 
#if I set gene id's to rownames, i cannot run the following command
#because there're duplicates
#after checking the duplicates, they are always two duplicates and one row has always zero overall
# THIS STEP WILL NOT BE NECESSARY SINCE I'M USING SYMBOLS GENERATED FROM THE GTF FILES

###drop all unneccessary cols
human_MOF <- human_MOF %>% dplyr::select(Geneid, contains("SRR"))
dim(human_MOF) #67062    22 
# 21 human samples - day 0, 2, 3, 5, 10, 15, 25 - 3 technical replicates
gorilla_MOF <- gorilla_MOF %>% dplyr::select(Geneid, contains("SRR"))
dim(gorilla_MOF) #67062    22 
# 21 gorilla samples - day 0, 2, 3, 5, 10, 15, 25 - 3 technical replicates

###how many genes are expressed in each species separately?
prefilter_h <- rowSums(human_MOF[2:22]) >= 10
cts_h_filt <- human_MOF[prefilter_h,] 
dim(cts_h_filt) #37001    22 - 30061 rows are removed
prefilter_g <- rowSums(gorilla_MOF[2:22]) >= 10
cts_g_filt <- gorilla_MOF[prefilter_g,] 
dim(cts_g_filt) #34826    22 - 32236 rows are removed

### merge human and non-human primate data together
cts <- cbind(human_MOF, gorilla_MOF %>% dplyr::select(!Geneid))
dim(cts) #67062 43

### import the table containing gene ids and symbols 
###which are extracted from "gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz"
ID_symbol <- read.table("../../Volumes/Extreme SSD/PhD/Drylab/Lancaster/data/ids_names_extractedFromGTF.txt")
colnames(ID_symbol) <- c("Geneid", "Genesymbol")

#In case of removing ensembl id version
#sub("\\.\\d+", "", rownames(cts_norm))

### add a new column of gene symbol using the imported table above
cts <- cbind(ID_symbol$Genesymbol, cts)
colnames(cts)[1] <- "Genesymbol"

sum(duplicated(cts$Genesymbol)) #4433 duplications
sum(duplicated(cts$Geneid)) #0 duplication

### remove rows having "_PAR_Y" in gene id
cts_mod <- cts[!grepl("_PAR_Y", cts$Geneid),] #45 gene ids contain "_PAR_Y"
dim(cts_mod) #67017    44

sum(duplicated(cts_mod$Genesymbol)) #4388 duplications
sum(duplicated(cts_mod$Geneid))

### pre-filter the data 
prefilter <- rowSums(cts_mod[3:44]) >= 10
cts_filt <- cts_mod[prefilter,] 
dim(cts_filt) #41204    44 - 25813 rows are removed

sum(duplicated(cts_filt$Genesymbol)) #28 duplications
sum(duplicated(cts_filt$Geneid))

###Preparing sample information table
x <- c(0,2,3,5,10,15,25)
h <- stri_paste("H9_", "d", x)
g <- stri_paste("G1_", "d", x)
Samples_h <- c()
Samples_g <- c()
for (i in 1:3){
  h_rep <- stri_paste(h, "_t", i)
  Samples_h <- c(Samples_h, h_rep)
  g_rep <- stri_paste(g, "_t", i)
  Samples_g <- c(Samples_g, g_rep)
  Samples <- c(Samples_h, Samples_g)
}

Samples <- c(Samples_h, Samples_g)

tp <- c(rep(stri_paste("day",x), 3))
sp <- c(rep("human", 21), rep("gorilla", 21))
coldata <- data.frame(tp, sp)
rownames(coldata) <- Samples

#Don't forget to change the colnames of count matrix
colnames(cts_filt)[3:44] <- Samples

### import the table containing gene ids and symbols 
###which are extracted from "gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz"
ID_symbol <- read.table(file = "data/ids_names_extractedFromGTF.txt")
colnames(ID_symbol) <- c("Geneid", "Genesymbol")

###import the table containing lncRNAs
lncRNA <- read.table(file = "../../Volumes/Extreme SSD/PhD/Drylab/Lancaster/data/GRFs_in_GRCh38_complete_list.txt", header = T)
lncRNA <- lncRNA %>% filter(Gene_type == "lncRNA")

#############################################################################
#DESeq2
#############################################################################
###data prep for DESeq2
cts_deseq <- cts_filt
rownames(cts_deseq) <- cts_deseq$Geneid
cts_deseq <- round(cts_deseq[,3:44])

write.table(cts_filt, file = "/Volumes/Extreme SSD/PhD/Drylab/Lancaster/output/preliminary/Lancaster_ReadCountTable_DESeq2_prefiltered.tsv", 
            sep = "\t", row.names = T, col.names = NA, quote = F)

#constructing a DESeqDataSet with multi-factor design
dds <- DESeqDataSetFromMatrix(countData = cts_deseq,
                              colData = coldata,
                              design = ~ tp + sp + tp:sp)
dds #dim: 41204 42 

#Setting the control samples
dds$tp <- factor(dds$tp, levels = c(stri_paste("day",x)))
dds$sp <- factor(dds$sp, levels = c("human", "gorilla"))

#DE analysis with full model 
dds <- DESeq(dds)

#DE analysis with reduced model (+ time-series)
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)
res_lrt
summary(res_lrt)
resultsNames(dds_lrt)

# Subset the LRT results to return genes with padj < 0.05
sig_res_LRT <- res_lrt %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)

# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
  pull(gene)
length(sigLRT_genes) #28174

# Obtain rlog values for those significant genes
cts_norm <- counts(dds, normalized=TRUE)
cts_cluster <- cts_norm[sig_res_LRT$gene, ]

#############################################################################
###run this part on the server
#############################################################################
saveRDS(cts_cluster, file = here("output/preliminary/cts_cluster.rds"))
saveRDS(coldata, file = here("output/preliminary/coldata.rds"))
#clusters <- degPatterns(cts_cluster, metadata = coldata, time = "tp", col=NULL)
#############################################################################
res <- results(dds)

cts_cluster <- readRDS("output/preliminary/cts_cluster.rds")
coldata <- readRDS("output/preliminary/coldata.rds")

saveRDS(cts_norm, file = here("output/preliminary/cts_norm_FINAL.rds"))
