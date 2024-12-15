#############################################################################
#DESeq2
#############################################################################
###data prep for DESeq2
cts_deseq <- cts_filt[,-11]
cts_deseq <- round(cts_deseq)

#Don't forget to change the colnames of count matrix
colnames(cts_deseq) <- Samples
cts_deseq[which(rownames(cts_deseq)  == "ENSG00000169554.23"),] 

#constructing a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = cts_deseq,
                              colData = coldata,
                              design = ~ sp)
dds #dim: 45492 10 
#constructing a DESeqDataSet with multi-factor design
#dds.n <- DESeqDataSetFromMatrix(countData = cts_deseq,
                                colData = coldata,
                                design = ~ sp + rep)
#dds.n #dim: 45492 10

#Setting the control samples
dds$sp <- factor(dds$sp, levels = c("human", "chimp"))
dds$rep <- factor(dds$rep, levels = c("rep1", "rep2"))
dds$Brep <- factor(dds$rep, levels = c("H1", "H2", "H3", "C1", "C2"))

#DE analysis with full model 
dds <- DESeq(dds)
resultsNames(dds) #"sp_chimp_vs_human"

#dds.n <- DESeq(dds.n)
#resultsNames(dds.n)

#PCA
vsd <- vst(dds, blind=FALSE)

#png("output/DEseq2/2024_09_30_PCA.png", width = 350, height = 350)

plotPCA(vsd, intgroup=c("sp")) + #PCA with species
  guides(color=guide_legend("species")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(plot.title = element_text(face="bold", size = 16)) +
  geom_point(size = 2) 

 #dev.off() 

#############################################################################
#ZEB2 exp
plotCounts(dds, gene= "ENSG00000169554.23", intgroup="sp")
#noramized read count table
cts_norm <- counts(dds, normalized=TRUE)
#cts_norm.n <- counts(dds.n, normalized=TRUE)

#ZEB2 exp per species
#png("output/DEseq2/2024_08_15_ZEB2exp_perSpecies.png", width = 350, height = 350)

cts_norm %>% as.data.frame() %>% filter(rownames(cts_norm) == "ENSG00000169554.23") %>%
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>%
  #mutate(rep = rep) %>%
  mutate(species = sp) %>%
  #mutate(rep = factor(rep, levels = c("rep1", "rep2"))) %>%
  mutate(species = factor(sp, levels = c("human", "chimp"))) %>%
  ggplot(aes(x = species, y = count, color = species)) +
  geom_boxplot() + stat_summary(fun.y = mean, geom = "line") +
  stat_compare_means(label = "p.signif", method = "wilcox.test", hide.ns = F) +
  scale_color_manual(values=c("#8dd35fff", "#ff7fb0ff")) +
  scale_y_log10() + theme_bw() + labs(title="Counts of ZEB2 in each species", x="") +
  ylab("normalized count (log)")

#dev.off()

#ZEB2 exp per biological replicates & species
#png("output/DEseq2/2024_08_15_ZEB2exp_perSpecies_perBrep.png", width = 350, height = 350)

cts_norm %>% as.data.frame() %>% filter(rownames(cts_norm) == "ENSG00000169554.23") %>%
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>%
  mutate(Brep = Brep) %>%
  mutate(species = sp) %>%
  mutate(Brep = factor(Brep, levels = c("H1", "H2", "H3", "C1", "C2"))) %>%
  mutate(species = factor(sp, levels = c("human", "chimp"))) %>%
  ggplot(aes(x = Brep, y = count, color = species, group = Brep)) +
  geom_boxplot() + stat_summary(fun.y = mean, geom = "line") +
  scale_color_manual(values=c("#8dd35fff", "#ff7fb0ff")) +
  scale_y_log10() + theme_bw() + labs(title="Counts of ZEB2 in each species", x="") +
  ylab("normalized count (log)")

#dev.off()
#############################################################################
#Log fold change shrinkage
resLFC <- lfcShrink(dds, coef="sp_chimp_vs_human", type="apeglm")
resLFC
plotMA(resLFC)
#ZEB2
resLFC[which(rownames(resLFC) == "ENSG00000169554.23"),]

#Final result
res0.01 <- resLFC[which(resLFC$padj < 0.01), ] #4176 DEGs
summary(res0.01)
saveRDS(res0.01, file = "preliminary/DESeq2_res_adjP_filtered.rds")

path_DEseq2 <- "output/DEseq2/"
write.table(ID_symbol %>% filter(Geneid %in% rownames(res0.01)[which(res0.01$log2FoldChange > 0)]), 
            file = paste(path_DEseq2,
                            'DEGs_posinCcomparedtoH_n', 
                            length(which(res0.01$log2FoldChange > 0)),
                            '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")
            
write.table(ID_symbol %>% filter(Geneid %in% rownames(res0.01)[which(res0.01$log2FoldChange < 0)]), 
            file = paste(path_DEseq2,
                            'DEGs_neginCcomparedtoH_n', 
                            length(which(res0.01$log2FoldChange < 0)),
                            '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

