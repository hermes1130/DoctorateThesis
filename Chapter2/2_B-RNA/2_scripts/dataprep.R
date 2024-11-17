library(DESeq2)
library(dplyr)
library(ggplot2)
library(here)
library(tidyr)
library(biomaRt)
library(psych)
library(stringi)
library(stringr)

#############################################################################
#Data prep
#############################################################################
### read counts preparation
human_rc <- read.table(here("../data/readcounts/GENCODE_incUnknown_genecounts_HSref"), header = T)
dim(human_rc) #67062 45

chimp_rc <- read.table(here("../data/readcounts/GENCODE_incUnknown_genecounts_Chimp_pantro6_crsmapped"), header = T)
dim(chimp_rc) #67062 18

orang_rc <- read.table(here("../data/readcounts/GENCODE_incUnknown_genecounts_Orang_ponabe3_crsmapped"), header = T)
dim(orang_rc) #67062 18

###drop all unneccessary cols
human_rc <- human_rc %>% dplyr::select(Geneid, contains("GM"))
dim(human_rc) #67062 10

chimp_rc <- chimp_rc %>% dplyr::select(Geneid, contains("Chimp")) %>% dplyr::select(!contains("fibro"))
dim(chimp_rc) #67062 10

orang_rc <- orang_rc %>% dplyr::select(Geneid, contains("Orang")) %>% dplyr::select(!contains("fibro"))
dim(orang_rc) #67062 10


### merge human and non-human primate data together
cts <- cbind(human_rc, chimp_rc %>% dplyr::select(!Geneid), orang_rc %>% dplyr::select(!Geneid))
dim(cts) #67062 28

### import the table containing gene ids and symbols 
###which are extracted from "gencode.v33.chr_patch_hapl_scaff.annotation.gtf.gz"
ID_symbol <- read.table(here("../data/ids_names_extractedFromGTF.txt"))
colnames(ID_symbol) <- c("Geneid", "Genesymbol")

### add a new column of gene symbol using the imported table above
cts <- cbind(ID_symbol$Genesymbol, cts)
colnames(cts)[1] <- "Genesymbol"

sum(duplicated(cts$Genesymbol)) #4433 duplications
sum(duplicated(cts$Geneid)) #0 duplication

### remove rows having "_PAR_Y" in gene id
cts_mod <- cts[!grepl("_PAR_Y", cts$Geneid),] #45 gene ids contain "_PAR_Y"
dim(cts_mod) #67017 29

sum(duplicated(cts_mod$Genesymbol)) #4388 duplications
sum(duplicated(cts_mod$Geneid)) #0

### pre-filter the data 
prefilter <- rowSums(cts_mod[3:29]) >= 10
cts_filt <- cts_mod[prefilter,] 
dim(cts_filt) #33896 29 - 33121 rows are removed

sum(duplicated(cts_filt$Genesymbol)) #13 duplications
sum(duplicated(cts_filt$Geneid))

### meta data prep
samples <- c("Human_R1_zeb2_2", "Human_R1_zeb2_3", "Human_R1_zeb2_neg",
"Human_R2_zeb2_2", "Human_R2_zeb2_3", "Human_R2_zeb2_neg",
"Human_R3_zeb2_2", "Human_R3_zeb2_3", "Human_R3_zeb2_neg",
"Chimp_Judith_Zeb2_2", "Chimp_Judith_Zeb2_3", "Chimp_Judith_Zeb2_neg",
"Chimp_Leo_Zeb2_2", "Chimp_Leo_Zeb2_3", "Chimp_Leo_Zeb2_neg",
"Chimp_Maryke_Zeb2_2", "Chimp_Maryke_Zeb2_3", "Chimp_Maryke_Zeb2_neg",
"Orang_Guchi_Zeb2_2", "Orang_Guchi_Zeb2_3", "Orang_Guchi_Zeb2_neg",
"Orang_Jaqo_Zeb2_2", "Orang_Jaqo_Zeb2_3", "Orang_Jaqo_Zeb2_neg",
"Orang_JingJing_Zeb2_2", "Orang_JingJing_Zeb2_3", "Orang_JingJing_Zeb2_neg")
#condition <- c(rep(c("KD1", "KD2", "control"), 9)) #separate KD
condition <- c(rep(c(rep("KD",2), "control"), 9)) #merged KD
species <- c(rep("Human", 9), rep("Chimpanzee", 9), rep("Orangutan", 9))

coldata <- data.frame(condition, species) 
rownames(coldata) <- samples

colnames(cts_filt)[3:29] <- samples

BG <- cts_filt[, 1:2]
write.table(BG, file = "output/preliminary/BackgroundGenes_forGO_allexpressedgenes_afterFilter.tsv", 
            sep = "\t", row.names = T, col.names = NA, quote = F)
