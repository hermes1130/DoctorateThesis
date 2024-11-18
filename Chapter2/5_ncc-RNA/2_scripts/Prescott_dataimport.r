library(here)
library(dplyr)
library(stringi)

library(DESeq2)
library(tibble)
library(DEGreport)

library(stringr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(plyr)
library(webr)
#the aim of this project: getting a list of DE genes between human and chimp NCCs


#############################################################################
#Data prep
#############################################################################
#read counts preparation
human_MOF <- read.table("data/readcounts/GENCODE_incUnknown_genecounts_human_hg38_MOF_Prescott", header = T)
dim(human_MOF) #70611    12
chimp_MOF <- read.table("data/readcounts/GENCODE_incUnknown_genecounts_chimp_hg38_MOF_Prescott", header = T)
dim(chimp_MOF) #70611    10

###version number drop 
#if I set gene id's to rownames, i cannot run the following command
#because there're duplicates
#after checking the duplicates, they are always two duplicates and one row has always zero overall
# THIS STEP WILL NOT BE NECESSARY SINCE I'M USING SYMBOLS GENERATED FROM THE GTF FILES

###drop all unneccessary cols
human_MOF <- human_MOF %>% dplyr::select(Geneid, contains("SRR"))
dim(human_MOF) #70611     7 
# 6 human samples - 3 biological replicates - 2 technical replicates
chimp_MOF <- chimp_MOF %>% dplyr::select(Geneid, contains("SRR"))
dim(chimp_MOF) #70611     5
# 4 chimp samples - 2 biological replicates - 2 technical replicates
#merge two species
cts <- cbind(human_MOF, chimp_MOF[,-1])

### import the table containing gene ids and symbols 
###which are extracted from "gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz"
ID_symbol <- read.table("data/ids_names_extractedFromGTF_v46.txt")
colnames(ID_symbol) <- c("Geneid", "Genesymbol")
ID_symbol[which(ID_symbol$Genesymbol  == "ZEB2"),] #ENSG00000169554.23

#In case of removing ensembl id version
#sub("\\.\\d+", "", rownames(cts_norm))

### add a new column of gene symbol using the imported table above
cts <- merge(cts,ID_symbol, by  = "Geneid") 

sum(duplicated(cts$Genesymbol)) #7244 duplications
sum(duplicated(cts$Geneid)) #0 duplication
rownames(cts) <- cts$Geneid
cts <- cts[,-1]
dim(cts) #70611    11

### pre-filter the data 
prefilter <- rowSums(cts[1:10]) >= 10
cts_filt <- cts[prefilter,] 
dim(cts_filt) #45492    11 - 25119 rows are removed

sum(duplicated(cts_filt$Genesymbol)) #618 duplications
sum(duplicated(cts_filt$Geneid))

write.table(cts_filt, file = "preliminary/Prescott_ReadCountTable_DESeq2_prefiltered.tsv", 
            sep = "\t", row.names = T, col.names = NA, quote = F)

###Preparing sample information table
hx <- c(1,2,3)
cx <- c(1,2)
h <- stri_paste("H", rep(hx,each = 2))
c <- stri_paste("C", rep(cx,each = 2))
Samples_h <- stri_paste(h, "_rep", c(1,2))
Samples_c <- stri_paste(c, "_rep", c(1,2))

Samples <- c(Samples_h, Samples_c)

sp <- c(rep("human", 6), rep("chimp", 4))
rep <- c(rep(stri_paste("rep",c(1,2)), 5))
Brep <- c(h, c)
coldata <- data.frame(sp, rep, Brep)
rownames(coldata) <- Samples

coldata$sp <- factor(coldata$sp)
coldata$rep <- factor(coldata$rep)
coldata$Brep <- factor(coldata$Brep)

#subject.n is the the individuals nested within a group, in my case species
#coldata$subject.n <- rep(c("1","2", "3", "1", "2"), each = 2)

###import the table containing lncRNAs
lncRNA <- read.table("data/GRFs_in_GRCh38_complete_list.txt", header = T)
lncRNA <- lncRNA %>% filter(Gene_type == "lncRNA")


