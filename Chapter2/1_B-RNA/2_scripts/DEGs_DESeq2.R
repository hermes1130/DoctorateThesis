library(DESeq2)
library(apeglm)
library(dplyr)
library(tibble)
library(ggbreak)
library(ggpubr)
library(pheatmap)
library(EnhancedVolcano)

library(plyr)

#############################################################################
#Data prep for DESeq2
#############################################################################
cts_deseq <- cts_filt
rownames(cts_deseq) <- cts_deseq$Geneid
cts_deseq <- cts_deseq[,3:29]
dim(cts_deseq) #33896    27

cts_deseq %>% as.data.frame() %>% filter(rownames(cts_deseq) == "ENSG00000169554.21") #ZEB2 is there

write.table(cts_deseq, file = "output/preliminary/ReadCountTable_DESeq2_prefiltered_allSpecies.tsv", 
            sep = "\t", row.names = T, col.names = NA, quote = F)

#############################################################################
#Run DESeq2
#############################################################################
### run dds
###this command is the standard argument for multiple design without pairing
dds <- DESeqDataSetFromMatrix(countData = cts_deseq,
                              colData = coldata,
                              design = ~ species + condition + species:condition)

###this command is for group-specific condition effects, individuals nested within groups
###this means each individuals are defined in coldata
#subject is the individual names
coldata$subject <- rep(c("R1","R2","R3","Judith", "Leo", "Maryke", "Guchi", "Jaqo", "JingJing"), each = 3)
#subject.n is the the individuals nested within a group, in my case species
coldata$subject.n <- rep(c("1","2", "3"),3, each = 3)

dds <- DESeqDataSetFromMatrix(countData = cts_deseq,
                              colData = coldata,
                              design = ~ species + species:subject.n + species:condition)


### note on factor levels
dds$condition #Levels: control KD
dds$species <- factor(dds$species, levels = c("Human","Chimpanzee","Orangutan"))
dds$subject # no levels defined
dds$subject.n #Levels: 1 2 3
  
### Differential expression analysis
dds <- DESeq(dds)
results(dds)
resultsNames(dds)
cts_norm <-counts(dds, normalized=TRUE)


#############################################################################
###main effect (KD effect in each species)
#############################################################################
#Several tries 
#1: without nested individuals within group & results() by subsetting by species: NO
#2: without nested individuals within group & results() true results, but doesn't count the individual differences in: NO
#3: with nested individuals within group & results() : YES

#############################################################################
###1: The very first try which doesn't really correspond to the read counts
res_cond_human <- results(dds[,dds$species == "Human"], name = "condition_KD_vs_control")
res_cond_chimp <- results(dds[,dds$species == "Chimpanzee"], name = "condition_KD_vs_control")
res_cond_orang <- results(dds[,dds$species == "Orangutan"], name = "condition_KD_vs_control")
summary(res_cond_human)
summary(res_cond_chimp)
summary(res_cond_orang)
res_cond_humanSig <- subset(res_cond_human, padj < 0.05)
res_cond_chimpSig <- subset(res_cond_chimp, padj < 0.05)
res_cond_orangSig <- subset(res_cond_orang, padj < 0.05)
sum(res_cond_humanSig$log2FoldChange > 0) #46
sum(res_cond_humanSig$log2FoldChange < 0) #43
sum(res_cond_chimpSig$log2FoldChange > 0) #15
sum(res_cond_chimpSig$log2FoldChange < 0) #12
sum(res_cond_orangSig$log2FoldChange > 0) #15
sum(res_cond_orangSig$log2FoldChange < 0) #12
#############################################################################
###2: The second try of the KD effect in each species - results are different from above
###Still doesn't correspond to the read counts
res_cond_human_test <- results(dds, contrast=c("condition","KD","control"))
summary(res_cond_human_test)
summary(subset(res_cond_human_test, padj < 0.05))
res_cond_humanSig_test <- subset(res_cond_human_test, padj < 0.05)

res_cond_chimp_test <- results(dds, list(c("condition_KD_vs_control", "speciesChimpanzee.conditionKD")))
summary(res_cond_chimp_test)
summary(subset(res_cond_chimp_test, padj < 0.05))
res_cond_chimpSig_test <- subset(res_cond_chimp_test, padj < 0.05)

res_cond_orang_test <- results(dds, list(c("condition_KD_vs_control", "speciesOrangutan.conditionKD")))
summary(res_cond_orang_test)
summary(subset(res_cond_orang_test, padj < 0.05))
res_cond_orangSig_test <- subset(res_cond_orang_test, padj < 0.05)
#############################################################################
###3: The third try with nested individuals within species

Main_H <- results(dds, name = "speciesHuman.conditionKD")
Main_H_LFC <- lfcShrink(dds, coef="speciesHuman.conditionKD", type="apeglm")
Main_C <- results(dds, name = "speciesChimpanzee.conditionKD")
Main_C_LFC <- lfcShrink(dds, coef="speciesChimpanzee.conditionKD", type="apeglm")
Main_O <- results(dds, name = "speciesOrangutan.conditionKD")
Main_O_LFC <- lfcShrink(dds, coef="speciesOrangutan.conditionKD", type="apeglm")
#schrinking doesn't look very convincing

Main_H_sig <- subset(Main_H, padj < 0.05)
Main_C_sig <- subset(Main_C, padj < 0.05)
Main_O_sig <- subset(Main_O, padj < 0.05)

NrDEs <- c(sum(Main_H_sig$log2FoldChange > 0),sum(Main_H_sig$log2FoldChange < 0), 
           sum(Main_C_sig$log2FoldChange > 0),sum(Main_C_sig$log2FoldChange < 0),
           sum(Main_O_sig$log2FoldChange > 0),sum(Main_O_sig$log2FoldChange < 0))

NrDEs <- NrDEs %>% as_tibble() %>% 
  add_column(species = rep(c("Human", "Chimpanzee", "Orangutan"), each =2 )) %>%
  add_column(reg = rep(c("UP", "DW"), 3)) 
NrDEs$value <- as.numeric(NrDEs$value)
NrDEs$species <- factor(NrDEs$species, levels = c("Human", "Chimpanzee", "Orangutan"))
NrDEs$reg <- factor(NrDEs$reg, levels = c("UP", "DW"))
NrDEs_sort <- arrange(NrDEs, species, reg)
NrDEs_sum <- ddply(NrDEs_sort, c("species"),
                     transform, label_ypos=cumsum(value))
NrDEs_sum$label_ypos <- c(97, 78,76, 68,626, 560)

#Number of DE genes
ggplot(data=NrDEs_sum, aes(x=species, y=value, fill=reg)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=value), vjust=1.1, color="black", size=3.0)+
  scale_y_continuous(limits=c(0, 650), breaks=c(0, 100, 600))+
  scale_y_break(c(150, 550)) +
  theme(axis.title = element_blank()) +
  labs(title="Number of DE genes (KD vs control)", x="", y="Number of DE genes") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2")

#if the group_by doesn't work, unload "plyr" 
NrDEs_per <- NrDEs_sum %>% 
  group_by(species) %>% 
  mutate(percentage = round(value/sum(value)*100))
NrDEs_per <- ddply(NrDEs_per, .(species),
                     transform, label_ypos_per = cumsum(percentage))
NrDEs_per$label_ypos_per <- c(88, 40, 92.5, 44.5,92.5,44.5)
  
#Percentage of up/down regulated genes
ggplot(data=NrDEs_per, aes(x=species, y=percentage, fill=reg)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos_per, label=percentage), vjust=0, color="black", size=4)+
  theme(axis.title = element_blank()) +
  labs(title="Composition of DE genes direction \n(KD vs control)", x="", y="Percentage (%)") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2")

#overlapping genes between species
intersect(rownames(Main_H_sig), rownames(Main_C_sig))
plotCounts(dds, gene="ENSG00000255319.5", intgroup=c("species", "condition"))
plotCounts(dds, gene="ENSG00000222881.1", intgroup=c("species", "condition"))
plotCounts(dds, gene="ENSG00000267215.2", intgroup=c("species", "condition"))
intersect(rownames(Main_H_sig), rownames(Main_O_sig)) #zero
intersect(rownames(Main_C_sig), rownames(Main_O_sig)) #zero

###boxplot of normalized read counts
cts_transf <- cts_norm %>% as.data.frame() %>% 
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>%
  mutate(condition = rep(coldata$condition, nrow(cts_norm))) %>%
  mutate(species = rep(coldata$species, nrow(cts_norm))) %>%
  mutate(gene = rep(rownames(cts_norm), each = 27))
cts_transf$species <- factor(cts_transf$species, levels=c("Human", "Chimpanzee", "Orangutan"))

cts_transf %>% filter(gene == "ENSG00000169554.21") %>%
  ggplot(aes(x = condition,y = count, fill = species)) +
  geom_boxplot() + 
  facet_wrap( ~ species) + theme_bw() +
  scale_fill_manual(values=c("#8dd35fff", "#ff80b2ff", "#ff9955ff")) +
  theme(axis.title = element_blank()) +
  labs(title="Normalized counts - ENSG00000169554.21", x="")
#############################################################################
###heatmap with all sig DE genes showing up/down regulation
#############################################################################
#this method is what suggested in deseq2 manual
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

#all DE genes detected in each species to be impacted significantly by ZEB2 KD
select <- which(rownames(dds) %in% rownames(mat))
label_deseq2heat <- as.data.frame(colData(dds)[,c("condition","species")])

pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=label_deseq2heat,
         cutree_rows = 8)
pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=label_deseq2heat,
         cutree_cols = 4, cutree_rows = 8)

#28 h-specific genes
select_hspc <- which(rownames(dds) %in% rownames(LFC_hsp))
label_deseq2heat_hspc <- as.data.frame(colData(dds)[,c("condition","species")])
pheatmap(assay(ntd)[select_hspc,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=label_deseq2heat_hspc)
#very low exp lvl

#############################################################################
#this method is based on the log2FC from each result table created by me
#############################################################################
#ALL genes
LFC_1 <- cbind(Main_H_sig[,2], 
               Main_C[rownames(Main_C) %in% rownames(Main_H_sig),2],
               Main_O[rownames(Main_O) %in% rownames(Main_H_sig),2])
rownames(LFC_1) <- rownames(Main_H_sig)
LFC_1.z <- t(apply(LFC_1, 1, scale))
colnames(LFC_1.z) <- c("Human", "Chimpanzee", "Orangutan")
pheatmap(LFC_1.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE)

LFC_2 <- cbind(Main_H[rownames(Main_H) %in% rownames(Main_C_sig),2],
               Main_C_sig[,2], 
               Main_O[rownames(Main_O) %in% rownames(Main_C_sig),2])
rownames(LFC_2) <- rownames(Main_C_sig)
LFC_2.z <- t(apply(LFC_2, 1, scale))
colnames(LFC_2.z) <- c("Human", "Chimpanzee", "Orangutan")
pheatmap(LFC_2.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE)

LFC_3 <- cbind(Main_H[rownames(Main_H) %in% rownames(Main_O_sig),2],
               Main_C[rownames(Main_C) %in% rownames(Main_O_sig),2],
               Main_O_sig[,2])
rownames(LFC_3) <- rownames(Main_O_sig)
LFC_3.z <- t(apply(LFC_3, 1, scale))
colnames(LFC_3.z) <- c("Human", "Chimpanzee", "Orangutan")
pheatmap(LFC_3.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE)

LFC_hsp <- cbind(Main_H[rownames(Main_H) %in% intersect(rownames(condDiff_HvsC_sig), rownames(condDiff_HvsO_sig)),2],
                 Main_C[rownames(Main_C) %in% intersect(rownames(condDiff_HvsC_sig), rownames(condDiff_HvsO_sig)),2],
                 Main_O[rownames(Main_O) %in% intersect(rownames(condDiff_HvsC_sig), rownames(condDiff_HvsO_sig)),2])
rownames(LFC_hsp) <- intersect(rownames(condDiff_HvsC_sig), rownames(condDiff_HvsO_sig))
LFC_hsp.z <- t(apply(LFC_hsp, 1, scale))
colnames(LFC_hsp.z) <- c("Human", "Chimpanzee", "Orangutan")
pheatmap(LFC_hsp.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE)
LFC_hsp.z_label <- LFC_hsp.z
rownames(LFC_hsp.z_label) <- ID_symbol %>% filter(Geneid %in% rownames(LFC_hsp.z_label)) %>% pull(Genesymbol)
pheatmap(LFC_hsp.z_label, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, color =rev(hcl.colors(30, palette = "GnBu")))

#15 protein-coding genes in 28 human-specific genes
pc15 <- rownames(LFC_hsp.z_label)[!str_detect(rownames(LFC_hsp.z_label), "^LINC|^BX|^AC|^AL|SEC13P1")]
LFC_15pc <- LFC_hsp.z_label[pc15,]
pheatmap(LFC_15pc, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE)


mat <- rbind(LFC_1, LFC_2, LFC_3)
mat.z <- t(apply(mat, 1, scale)) #scaling changes the direction of +/-
colnames(mat.z) <- c("Human", "Chimpanzee", "Orangutan")

#no gene names
label_whole <- ID_symbol %>% filter(Geneid %in% rownames(mat.z)) %>% pull(Genesymbol)
pheatmap(mat.z, cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=FALSE, labels_row = label_whole)

#three genes = three rows are duplicated
mat.z_label <- mat.z
nrow(mat.z_label)
mat.z_label_uniq <- mat.z_label[-which(duplicated(rownames(mat.z_label))),]
nrow(mat.z_label_uniq)
rownames(mat.z_label_uniq) <- ID_symbol %>% filter(Geneid %in% rownames(mat.z_label_uniq)) %>% pull(Genesymbol)
pheatmap(mat.z_label_uniq, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE)

#gene names with sig DEs from humans - wrong
label_h <- ID_symbol %>% filter(Geneid %in% rownames(Main_H_sig)) %>% pull(Genesymbol)
label_whole[!label_whole %in% label_h] <- "" 
pheatmap(mat.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, labels_row = label_whole)



#TOP10 from each 
Main_H_sigOrdered <- Main_H_sig[order(Main_H_sig$padj),]
Main_C_sigOrdered <- Main_C_sig[order(Main_C_sig$padj),]
Main_O_sigOrdered <- Main_O_sig[order(Main_O_sig$padj),]

#need to order the rownames of C and O by H!!!!
LFC_1_ordered <- cbind(Main_H_sigOrdered[1:10,2], 
               Main_C[match(rownames(Main_H_sigOrdered)[1:10], rownames(Main_C)),2],
               Main_O[match(rownames(Main_H_sigOrdered)[1:10], rownames(Main_O)),2])
rownames(LFC_1_ordered) <- rownames(Main_H_sigOrdered)[1:10]

LFC_2_ordered <- cbind(Main_H[match(rownames(Main_C_sigOrdered)[1:10], rownames(Main_H)),2],
               Main_C_sigOrdered[1:10,2], 
               Main_O[match(rownames(Main_C_sigOrdered)[1:10], rownames(Main_O)),2])
rownames(LFC_2_ordered) <- rownames(Main_C_sigOrdered)[1:10]

LFC_3_ordered <- cbind(Main_H[match(rownames(Main_O_sigOrdered)[1:10], rownames(Main_H)),2],
                       Main_C[match(rownames(Main_O_sigOrdered)[1:10], rownames(Main_C)),2],
               Main_O_sigOrdered[1:10,2])
rownames(LFC_3_ordered) <- rownames(Main_O_sigOrdered)[1:10]

mat_top10 <- rbind(LFC_1_ordered, LFC_2_ordered, LFC_3_ordered)
mat_mat_top10.z <- t(apply(mat_top10, 1, scale)) #scaling changes the direction of +/-
colnames(mat_mat_top10.z) <- c("Human", "Chimpanzee", "Orangutan")

pheatmap(mat_mat_top10.z, cluster_rows=TRUE, show_rownames=TRUE, show_colnames=TRUE,
         cluster_cols=FALSE, labels_row = ID_symbol %>% filter(Geneid %in% rownames(mat_mat_top10.z)) %>% pull(Genesymbol))

#############################################################################
#heatmap of 82 common peaks and 438 human-specific peaks
#############################################################################
CommonP <- read.table(file = "../../ZEB2_ChIPseq/ChIP-seq_results/data/Genes close to peaks - Overlap within 10k.tsv")
nrow(CommonP) #82
HscpP <- read.table(file = "../../ZEB2_ChIPseq/ChIP-seq_results/data/Genes close to peaks - Human specific targets within 10k.tsv")
nrow(HscpP) #437

CommonP_id <- ID_symbol %>% filter(Genesymbol %in% CommonP$V1) 
nrow(CommonP_id) #86
# last 4 gene symbols are duplicated since they have two diff ensembl IDs -.-

which(duplicated(ID_symbol %>% filter(Genesymbol %in% HscpP$V1) %>% pull(Genesymbol)))
#16 genes are duplicated
HscpP_id <- ID_symbol %>% filter(Genesymbol %in% HscpP$V1) 
nrow(HscpP_id) #453
#############################################################################
#LFC matrix - per species
#for 1) common peaks & 2) human-specific peaks
#1) common peaks
LFC_cp <- cbind(Main_H[rownames(Main_H) %in% CommonP_id$Geneid,2],
               Main_C[rownames(Main_C) %in% CommonP_id$Geneid,2], 
               Main_O[rownames(Main_O) %in% CommonP_id$Geneid,2])
nrow(LFC_cp) #again 82 genes
rownames(LFC_cp) <- ID_symbol %>% filter(Geneid %in% rownames(Main_H[rownames(Main_H) %in% CommonP_id$Geneid,])) %>% pull(Genesymbol)
LFC_cp.z <- t(apply(LFC_cp, 1, scale))
colnames(LFC_cp.z) <- c("Human", "Chimpanzee", "Orangutan")
pheatmap(LFC_cp.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, cutree_rows = 6)

#1) common peaks after filtering
filter_cp <- ID_symbol %>% filter(Geneid %in% rownames(heat_cp_filter)) %>% pull(Genesymbol)
LFC_cp_f.z <- LFC_cp.z %>% as.data.frame() %>% filter(rownames(LFC_cp.z) %in% filter_cp)
pheatmap(LFC_cp_f.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, cutree_rows = 4)

#2) human-specific peaks
LFC_hscp <- cbind(Main_H[rownames(Main_H) %in% HscpP_id$Geneid,2],
                Main_C[rownames(Main_C) %in% HscpP_id$Geneid,2], 
                Main_O[rownames(Main_O) %in% HscpP_id$Geneid,2])
nrow(LFC_hscp) # 247 genes
rownames(LFC_hscp) <- ID_symbol %>% filter(Geneid %in% rownames(Main_H[rownames(Main_H) %in% HscpP_id$Geneid,])) %>% pull(Genesymbol)
LFC_hscp.z <- t(apply(LFC_hscp, 1, scale))
colnames(LFC_hscp.z) <- c("Human", "Chimpanzee", "Orangutan")
pheatmap(LFC_hscp.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, cutree_rows = 6)

#2) human-specific peaks after filtering
filter_hscp <- ID_symbol %>% filter(Geneid %in% rownames(heat_hscp_filter)) %>% pull(Genesymbol)
LFC_hscp_f.z <- LFC_hscp.z %>% as.data.frame() %>% filter(rownames(LFC_hscp.z) %in% filter_hscp)
pheatmap(LFC_hscp_f.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, cutree_rows = 5)

#############################################################################
#LFC matrix - per individuals
#for 1) common peaks & 2) human-specific peaks
#1) common peaks n = 82
heat_cp_ind <- data.frame(matrix(, nrow=nrow(heat_cp), ncol=9))
rownames(heat_cp_ind) <- rownames(heat_cp)
colnames(heat_cp_ind) <- c("Human_R1","Human_R2", "Human_R3", 
                           "Chimp_Judith", "Chimp_Leo","Chimp_Maryke", 
                           "Orang_Guchi", "Orang_Jaqo","Orang_JingJing")
heat_cp_1 <- heat_cp+1
for(i in 1:ncol(heat_cp_ind)){
  heat_cp_ind[,i] <- log2(((heat_cp_1[,1+(3*(i-1))]/heat_cp_1[,3+(3*(i-1))]) + 
                             (heat_cp_1[,2+(3*(i-1))]/heat_cp_1[,3+(3*(i-1))]))/2)
}

rownames(heat_cp_ind) <- ID_symbol %>% filter(Geneid %in% rownames(heat_cp_ind)) %>% pull(Genesymbol)
myscale <- function(x) {
  x <- ifelse(x < 0 , -x[x<0]/min(x), x[x>=0]/max(x))
}
heat_cp_ind.z <- t(apply(heat_cp_ind, 1, myscale))

label_ind <- data.frame(species = c(rep("Human",3), rep("Chimp",3), rep("Orang",3)))
rownames(label_ind) <- colnames(heat_cp_ind.z)
pheatmap(heat_cp_ind.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=label_ind, cutree_rows = 4) 

#2) human-specific peaks n = 247
heat_hscp_ind <- data.frame(matrix(, nrow=nrow(heat_hscp), ncol=9))
rownames(heat_hscp_ind) <- rownames(heat_hscp)
colnames(heat_hscp_ind) <- c("Human_R1","Human_R2", "Human_R3", 
                           "Chimp_Judith", "Chimp_Leo","Chimp_Maryke", 
                           "Orang_Guchi", "Orang_Jaqo","Orang_JingJing")
heat_hscp_1 <- heat_hscp+1
for(i in 1:ncol(heat_hscp_ind)){
  heat_hscp_ind[,i] <- log2(((heat_hscp_1[,1+(3*(i-1))]/heat_hscp_1[,3+(3*(i-1))]) + 
                             (heat_hscp_1[,2+(3*(i-1))]/heat_hscp_1[,3+(3*(i-1))]))/2)
}

rownames(heat_hscp_ind) <- ID_symbol %>% filter(Geneid %in% rownames(heat_hscp_ind)) %>% pull(Genesymbol)
heat_hscp_ind.z <- t(apply(heat_hscp_ind, 1, myscale))

pheatmap(heat_hscp_ind.z, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=label_ind, cutree_rows = 3) 



#############################################################################
#normalized read count matrix for 1) common peaks & 2) human-specific peaks
label <- as.data.frame(colData(dds)[,c("condition","species")])

#1) common peaks
select_cp <- which(rownames(dds) %in% CommonP_id$Geneid)
length(select_cp) #82
heat_cp <- assay(ntd)[select_cp,]
heat_cp_filter <- heat_cp %>% as.data.frame() %>% filter_all(any_vars(. > 8)) #%>% nrow() #70

#heatmap of all 82 common peaks
pheatmap(assay(ntd)[select_cp,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=label, cutree_rows = 4)
#heat of 70 common peaks that are expressed above 5 in all samples
breaksList = seq(0, 20, by = 1)
pheatmap(heat_cp_filter, 
         cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=label, cutree_rows = 4,
         #labels_row = ID_symbol %>% filter(Geneid %in% rownames(heat_cp_filter)) %>% pull(Genesymbol),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)

#2) human-specific peaks
select_hscp <- which(rownames(dds) %in% HscpP_id$Geneid)
length(select_hscp) #247
heat_hscp <- assay(ntd)[select_hscp,]
heat_hscp_filter <- heat_hscp %>% as.data.frame() %>% filter_all(any_vars(. > 8)) #%>% nrow() #50

#heatmap of all 247 human-specific peaks
pheatmap(assay(ntd)[select_hscp,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=label, cutree_rows = 4)
#heatmap of 50 human-specific peaks that are expressed above 8 in any samples
pheatmap(heat_hscp_filter, 
         cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=label,cutree_rows = 4,
         #labels_row = ID_symbol %>% filter(Geneid %in% rownames(heat_hscp_filter)) %>% pull(Genesymbol),
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList)


#############################################################################
#volcano plot
#############################################################################
EnhancedVolcano(Main_H,
                lab = ID_symbol[ID_symbol$Geneid %in% rownames(Main_H),2],
                x = 'log2FoldChange',
                y = 'padj',
                title = 'KD effect in HS',
                pCutoff = 10e-2/2,
                FCcutoff = 0,
                col=c('#b5deff', '#b5deff', '#b5deff', 'red3'),
                labSize = 2.0,
                pointSize = 0.5,
                axisLabSize = 10.0)

EnhancedVolcano(Main_C,
                lab = ID_symbol[ID_symbol$Geneid %in% rownames(Main_C),2],
                x = 'log2FoldChange',
                y = 'padj',
                title = 'KD effect in PT',
                pCutoff = 10e-2/2,
                FCcutoff = 0,
                col=c('#b5deff', '#b5deff', '#b5deff', 'red3'),
                labSize = 2.0,
                pointSize = 0.5,
                axisLabSize = 10.0)

EnhancedVolcano(Main_O,
                lab = ID_symbol[ID_symbol$Geneid %in% rownames(Main_O),2],
                x = 'log2FoldChange',
                y = 'padj',
                title = 'KD effect in PA',
                pCutoff = 10e-2/2,
                FCcutoff = 0,
                col=c('#b5deff', '#b5deff', '#b5deff', 'red3'),
                labSize = 2.0,
                pointSize = 0.5,
                axisLabSize = 10.0)



#############################################################################
#Output saving
#############################################################################
###1:
#up_H (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_humanSig[which(res_cond_humanSig$log2FoldChange > 0),])),
            file = here("output/results/KDeffect_UP_Human_n46.tsv"), 
            quote = F, row.names = F)
#down_H (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_humanSig[which(res_cond_humanSig$log2FoldChange < 0),])),
            file = here("output/results/KDeffect_DW_Human_n43.tsv"), 
            quote = F, row.names = F)
#up_C (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_chimpSig[which(res_cond_chimpSig$log2FoldChange > 0),])),
            file = here("output/results/KDeffect_UP_Chimp_n15.tsv"), 
            quote = F, row.names = F)
#down_C (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_chimpSig[which(res_cond_chimpSig$log2FoldChange < 0),])),
            file = here("output/results/KDeffect_DW_Chimp_n12.tsv"), 
            quote = F, row.names = F)
#up_O (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_orangSig[which(res_cond_orangSig$log2FoldChange > 0),])),
            file = here("output/results/KDeffect_UP_Orang_n15.tsv"), 
            quote = F, row.names = F)
#down_O (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_orangSig[which(res_cond_orangSig$log2FoldChange < 0),])),
            file = here("output/results/KDeffect_DW_Orang_n12.tsv"), 
            quote = F, row.names = F)


###2:output from the NEW APPROACH
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_humanSig_test[which(res_cond_humanSig_test$log2FoldChange > 0),])),
            file = here("../output/results/NEW_KDeffect_UP_Human_n15.tsv"), 
            quote = F, row.names = F)
#down_H (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_humanSig_test[which(res_cond_humanSig_test$log2FoldChange < 0),])),
            file = here("../output/results/NEW_KDeffect_DW_Human_n12.tsv"), 
            quote = F, row.names = F)
#up_C (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_chimpSig_test[which(res_cond_chimpSig_test$log2FoldChange > 0),])),
            file = here("../output/results/NEW_KDeffect_UP_Chimp_n5.tsv"), 
            quote = F, row.names = F)
#down_C (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_chimpSig_test[which(res_cond_chimpSig_test$log2FoldChange < 0),])),
            file = here("../output/results/NEW_KDeffect_DW_Chimp_n2.tsv"), 
            quote = F, row.names = F)
#up_O (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_orangSig_test[which(res_cond_orangSig_test$log2FoldChange > 0),])),
            file = here("../output/results/NEW_KDeffect_UP_Orang_n10.tsv"), 
            quote = F, row.names = F)
#down_O (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_cond_orangSig_test[which(res_cond_orangSig_test$log2FoldChange < 0),])),
            file = here("../output/results/NEW_KDeffect_DW_Orang_n441.tsv"), 
            quote = F, row.names = F)
###3:output from the nesting
#up_H (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(Main_H_sig[which(Main_H_sig$log2FoldChange > 0),])),
            file = here("/Volumes/Extreme SSD/PhD/Drylab/ZEB2_RNAseq/ZEB2_RNAseq_final/output/results/NESTED_KDeffect_UP_Human_n19.tsv"), 
            quote = F, row.names = F)
#down_H (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(Main_H_sig[which(Main_H_sig$log2FoldChange < 0),])),
            file = here("/Volumes/Extreme SSD/PhD/Drylab/ZEB2_RNAseq/ZEB2_RNAseq_final/output/results/NESTED_KDeffect_DW_Human_n78.tsv"), 
            quote = F, row.names = F)
#up_C (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(Main_C_sig[which(Main_C_sig$log2FoldChange > 0),])),
            file = here("/Volumes/Extreme SSD/PhD/Drylab/ZEB2_RNAseq/ZEB2_RNAseq_final/output/results/NESTED_KDeffect_UP_Chimp_n8.tsv"), 
            quote = F, row.names = F)
#down_C (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(Main_C_sig[which(Main_C_sig$log2FoldChange < 0),])),
            file = here("/Volumes/Extreme SSD/PhD/Drylab/ZEB2_RNAseq/ZEB2_RNAseq_final/output/results/NESTED_KDeffect_DW_Chimp_n68.tsv"), 
            quote = F, row.names = F)
#up_O (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(Main_O_sig[which(Main_O_sig$log2FoldChange > 0),])),
            file = here("/Volumes/Extreme SSD/PhD/Drylab/ZEB2_RNAseq/ZEB2_RNAseq_final/output/results/NESTED_KDeffect_UP_Orang_n66.tsv"), 
            quote = F, row.names = F)
#down_O (KD effect)
write.table(ID_symbol %>% filter(Geneid %in% rownames(Main_O_sig[which(Main_O_sig$log2FoldChange < 0),])),
            file = here("/Volumes/Extreme SSD/PhD/Drylab/ZEB2_RNAseq/ZEB2_RNAseq_final/output/results/NESTED_KDeffect_DW_Orang_n560.tsv"), 
            quote = F, row.names = F)


#############################################################################
###1: The first try without nesting the individuals
###main effect (KD effect) + interaction term (main effect in NHP compared to human)
res_chimp <- results(dds, contrast =list("condition_KD_vs_control", "speciesChimpanzee.conditionKD"))
res_orang <- results(dds, contrast =list("condition_KD_vs_control", "speciesOrangutan.conditionKD"))
res_c_o <- results(dds, contrast =list("speciesChimpanzee.conditionKD", "speciesOrangutan.conditionKD"))
summary(res_chimp)
summary(res_orang)
summary(res_c_o)

res_chimpSig <- subset(res_chimp, padj < 0.05)
res_orangSig <- subset(res_orang, padj < 0.05)
res_c_oSig <- subset(res_c_o, padj < 0.05)
sum(res_chimpSig$log2FoldChange > 0) #15
sum(res_chimpSig$log2FoldChange < 0) #6
sum(res_orangSig$log2FoldChange > 0) #15
sum(res_orangSig$log2FoldChange < 0) #5
sum(res_c_oSig$log2FoldChange > 0) #8
sum(res_c_oSig$log2FoldChange < 0) #3

###2: The second try with nested individuals within species
###Some tests
res_HvsC_1 <- results(dds, name = "species_Chimpanzee_vs_Human") 
#Species diff without considering condition
res_HvsC_2 <- results(dds, contrast = list("speciesChimpanzee.conditionKD", "speciesHuman.conditionKD"))
#condition effect difference between two species
#THIS ONE!!!

subset(res_HvsC_1, padj < 0.05) #4610
subset(res_HvsC_2, padj < 0.05) #70

###DE analysis
condDiff_HvsC <- results(dds, contrast = list("speciesChimpanzee.conditionKD", "speciesHuman.conditionKD"))
condDiff_HvsO <- results(dds, contrast = list("speciesOrangutan.conditionKD", "speciesHuman.conditionKD"))
condDiff_CvsO <- results(dds, contrast = list("speciesChimpanzee.conditionKD", "speciesOrangutan.conditionKD"))

nrow(subset(condDiff_HvsC, padj < 0.1)) #90
nrow(subset(condDiff_HvsO, padj < 0.1)) #146
nrow(subset(condDiff_CvsO, padj < 0.1)) #99

condDiff_HvsC_sig <- subset(condDiff_HvsC, padj < 0.05) #70
condDiff_HvsO_sig <- subset(condDiff_HvsO, padj < 0.05) #101
condDiff_CvsO_sig <- subset(condDiff_CvsO, padj < 0.05) #85

condDiff <- c(sum(condDiff_HvsC_sig$log2FoldChange > 0),sum(condDiff_HvsC_sig$log2FoldChange < 0), 
           sum(condDiff_HvsO_sig$log2FoldChange > 0),sum(condDiff_HvsO_sig$log2FoldChange < 0),
           sum(condDiff_CvsO_sig$log2FoldChange > 0),sum(condDiff_CvsO_sig$log2FoldChange < 0))

condDiff <- condDiff %>% as_tibble() %>% 
  add_column(species = rep(c("HvsC", "HvsO", "CvsO"), each =2 )) %>%
  add_column(reg = rep(c("UP", "DW"), 3)) 
condDiff$value <- as.numeric(condDiff$value)
condDiff$species <- factor(condDiff$species, levels = c("HvsC", "HvsO", "CvsO"))
condDiff$reg <- factor(condDiff$reg, levels = c("UP", "DW"))
condDiff_sort <- arrange(condDiff, species, reg)
condDiff_sum <- ddply(condDiff_sort, c("species"),
                   transform, label_ypos=cumsum(value))
condDiff_sum$label_ypos <- c(70, 39,101, 45,85, 48)

#Number of DE genes
ggplot(data=condDiff_sum, aes(x=species, y=value, fill=reg)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=value), vjust=1.1, color="black", size=3.0)+
  scale_y_continuous(limits=c(0, 110), breaks=c(0, 50, 100))+
  theme(axis.title = element_blank()) +
  labs(title="Number of DE genes \n(Condition effect difference)", x="", y="Number of DE genes") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2")

#overlapping genes between species
geneList1 <- intersect(rownames(condDiff_HvsC_sig), rownames(condDiff_HvsO_sig)) #28
geneList2 <- intersect(rownames(condDiff_HvsC_sig), rownames(condDiff_CvsO_sig)) #14
geneList3 <- intersect(rownames(condDiff_HvsO_sig), rownames(condDiff_CvsO_sig)) #19

###boxplot of normalized read counts
cts_transf <- cts_norm %>% as.data.frame() %>% 
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>%
  mutate(condition = rep(coldata$condition, nrow(cts_norm))) %>%
  mutate(species = rep(c("H", "C", "O"), each = 9, nrow(cts_norm))) %>%
  mutate(gene = rep(rownames(cts_norm), each = 27))
cts_transf$species <- factor(cts_transf$species, levels=c("H", "C", "O"))

cts_transf %>% filter(gene == geneList1[1]) %>%
  ggplot(aes(x = condition,y = count, fill = species)) +
  geom_boxplot() + 
  facet_grid(gene ~ species, scales="free_y") + theme_bw() +  theme(legend.position="none") +
  scale_fill_manual(values=c("#8dd35fff", "#ff80b2ff", "#ff9955ff")) +
  theme(axis.title = element_blank()) +
  theme(axis.text=element_text(size=5))

plotList <- lapply(
  1:length(geneList1),
  function(genename) {
    # Need to assign the plot to a variable because 
    # you want to generate the plot AND save to file 
    x <- cts_transf %>% filter(gene == geneList1[genename]) %>%
      ggplot(aes(x = condition,y = count, fill = species)) +
      geom_boxplot() + 
      facet_grid(gene ~ species, scales="free_y") + theme_bw() +  theme(legend.position="none") +
      scale_fill_manual(values=c("#8dd35fff", "#ff80b2ff", "#ff9955ff")) +
      theme(axis.title = element_blank()) +
      theme(axis.text=element_text(size=5))
    # No need for the plot argument.  It defaults to the last plot created.
    ggsave(filename = paste0("plots/", genename, ".jpg"))
    # Return the plot just created
    x
  }
)

allplots <- ggarrange(plotlist=plotList,
                      ncol = 7, nrow = 4)

allplots

#volcano plots
EnhancedVolcano(condDiff_HvsC,
                lab = ID_symbol[ID_symbol$Geneid %in% rownames(condDiff_HvsC),2],
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'KD effect changes in HvsC',
                pCutoff = 10e-2/2,
                FCcutoff = 0,
                col=c('#b5deff', '#b5deff', '#b5deff', 'red3'),
                labSize = 2.0,
                pointSize = 0.5,
                axisLabSize = 10.0)

#############################################################################
#Output saving
#############################################################################
#1:
#up_HvsC (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_chimpSig[which(res_chimpSig$log2FoldChange > 0),])),
            file = here("output/results/KDeffect_SpeciesDiff_UP_HumanvsChimp_n15.tsv"), 
            quote = F, row.names = F)
#down_HvsC (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_chimpSig[which(res_chimpSig$log2FoldChange < 0),])),
            file = here("output/results/KDeffect_SpeciesDiff_DW_HumanvsChimp_n6.tsv"), 
            quote = F, row.names = F)
#up_HvsO (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_orangSig[which(res_orangSig$log2FoldChange > 0),])),
            file = here("output/results/KDeffect_SpeciesDiff_UP_HumanvsOrang_n15.tsv"), 
            quote = F, row.names = F)
#down_HvsO (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_orangSig[which(res_orangSig$log2FoldChange < 0),])),
            file = here("output/results/KDeffect_SpeciesDiff_DW_HumanvsOrang_n5.tsv"), 
            quote = F, row.names = F)
#up_CvsO (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_c_oSig[which(res_c_oSig$log2FoldChange > 0),])),
            file = here("output/results/KDeffect_SpeciesDiff_UP_ChimpvsOrang_n8.tsv"), 
            quote = F, row.names = F)
#down_CvsO (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(res_c_oSig[which(res_c_oSig$log2FoldChange < 0),])),
            file = here("output/results/KDeffect_SpeciesDiff_DW_ChimpvsOrang_n3.tsv"), 
            quote = F, row.names = F)

#2:
#up_HvsC (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(condDiff_HvsC_sig[which(condDiff_HvsC_sig$log2FoldChange > 0),])),
            file = here("output/results/NESTED_KDeffect_SpeciesDiff_UP_HumanvsChimp_n31.tsv"), 
            quote = F, row.names = F)
#down_HvsC (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(condDiff_HvsC_sig[which(condDiff_HvsC_sig$log2FoldChange < 0),])),
            file = here("output/results/NESTED_KDeffect_SpeciesDiff_DW_HumanvsChimp_n39.tsv"), 
            quote = F, row.names = F)
#up_HvsO (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(condDiff_HvsO_sig[which(condDiff_HvsO_sig$log2FoldChange > 0),])),
            file = here("output/results/NESTED_KDeffect_SpeciesDiff_UP_HumanvsOrang_n56.tsv"), 
            quote = F, row.names = F)
#down_HvsO (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(condDiff_HvsO_sig[which(condDiff_HvsO_sig$log2FoldChange < 0),])),
            file = here("output/results/NESTED_KDeffect_SpeciesDiff_DW_HumanvsOrang_n45.tsv"), 
            quote = F, row.names = F)
#up_CvsO (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(condDiff_CvsO_sig[which(condDiff_CvsO_sig$log2FoldChange > 0),])),
            file = here("output/results/NESTED_KDeffect_SpeciesDiff_UP_ChimpvsOrang_n37.tsv"), 
            quote = F, row.names = F)
#down_CvsO (KD effect + species diff)
write.table(ID_symbol %>% filter(Geneid %in% rownames(condDiff_CvsO_sig[which(condDiff_CvsO_sig$log2FoldChange < 0),])),
            file = here("output/results/NESTED_KDeffect_SpeciesDiff_DW_ChimpvsOrang_n48.tsv"), 
            quote = F, row.names = F)

#h-specific
write.table(ID_symbol %>% filter(Geneid %in% rownames(condDiff_HvsC_sig[which(rownames(condDiff_HvsC_sig) %in% rownames(condDiff_HvsO_sig)),])),
            file = here("/Volumes/Extreme SSD/PhD/Drylab/ZEB2_RNAseq/ZEB2_RNAseq_final/output/results/NESTED_Hspec_n28.tsv"), 
            quote = F, row.names = F)

#c-specific
write.table(ID_symbol %>% filter(Geneid %in% intersect(rownames(condDiff_HvsC_sig), rownames(condDiff_CvsO_sig))),
            file = "/Volumes/Extreme SSD/PhD/Drylab/ZEB2_RNAseq/ZEB2_RNAseq_final/output/results/2024_09_06_NESTED_Cspec_n14.tsv", 
            quote = F, row.names = F, col.names = T)

#o-specific
write.table(ID_symbol %>% filter(Geneid %in% intersect(rownames(condDiff_HvsO_sig), rownames(condDiff_CvsO_sig))),
            file = "/Volumes/Extreme SSD/PhD/Drylab/ZEB2_RNAseq/ZEB2_RNAseq_final/output/results/2024_09_06_NESTED_Ospec_n19.tsv", 
            quote = F, row.names = F, col.names = T)

#############################################################################
### PCA
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("species")) + theme_bw() +scale_color_manual(values=c("#8dd35fff", "#ff80b2ff", "#ff9955ff"))
plotPCA(vsd, intgroup=c("condition")) + theme_bw()
plotPCA(vsd, intgroup=c("condition", "species")) + theme_bw()
plotPCA(vsd, intgroup=c("condition", "species", "subject.n")) + theme_bw()

plotPCA(vsd[,vsd$species == "Human"], intgroup=c("condition")) + theme_bw()
plotPCA(vsd[,vsd$species == "Chimpanzee"], intgroup=c("condition")) + theme_bw()
plotPCA(vsd[,vsd$species == "Orangutan"], intgroup=c("condition")) + theme_bw()


