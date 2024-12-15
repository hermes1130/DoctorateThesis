library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")

#############################################################################
#GO - CompareCluster
#############################################################################
#First try
#I want to create a emapplot with pos. and neg. correlated genes at each time point
genelist_forGO <- data.frame()

for (i in 1:ncol(corr_h)){
  h_pos <- rownames(corr_h)[which(corr_h[,i] > 0)] %>% 
    as.data.frame() %>% 
    mutate(TP = stri_paste("d", x)[i]) %>%
    mutate(reg = "(+)") %>%
    mutate(sp = "Human") %>%
    mutate(symbol = ID_symbol %>% 
             filter(Geneid %in% rownames(corr_h)[which(corr_h[,i] > 0)]) %>% 
             pull(2))
  h_neg <- rownames(corr_h)[which(corr_h[,i] < 0)] %>% 
    as.data.frame() %>% 
    mutate(TP = stri_paste("d", x)[i]) %>%
    mutate(reg = "(-)") %>%
    mutate(sp = "Human") %>%
    mutate(symbol = ID_symbol %>% 
             filter(Geneid %in% rownames(corr_h)[which(corr_h[,i] < 0)]) %>% 
             pull(2))
  
  g_pos <- rownames(corr_g)[which(corr_g[,i] > 0)] %>% 
    as.data.frame() %>% 
    mutate(TP = stri_paste("d", x)[i]) %>%
    mutate(reg = "(+)")  %>%
    mutate(sp = "Gorilla") %>%
    mutate(symbol = ID_symbol %>% 
             filter(Geneid %in% rownames(corr_h)[which(corr_g[,i] > 0)]) %>% 
             pull(2))
  g_neg <- rownames(corr_g)[which(corr_g[,i] < 0)] %>% 
    as.data.frame() %>% 
    mutate(TP = stri_paste("d", x)[i]) %>%
    mutate(reg = "(-)") %>%
    mutate(sp = "Gorilla") %>%
    mutate(symbol = ID_symbol %>% 
             filter(Geneid %in% rownames(corr_h)[which(corr_g[,i] < 0)]) %>% 
             pull(2))
  
  genelist_forGO <- rbind(genelist_forGO, h_pos, h_neg, g_pos, g_neg)
}
colnames(genelist_forGO)[1] <- "gene"

#compareCluster with the data above
compare.formular_all <- compareCluster(symbol~TP+reg, 
                                       data =genelist_forGO,
                                       pvalueCutoff = 0.05,
                                       fun = 'enrichGO',
                                       OrgDb= 'org.Hs.eg.db', 
                                       keyType = 'SYMBOL', 
                                       ont = 'BP',
                                       universe = cts_filt$Genesymbol)

compare.formular_HvsG <- compareCluster(symbol~sp+TP+reg, 
                                       data =genelist_forGO,
                                       pvalueCutoff = 0.05,
                                       fun = 'enrichGO',
                                       OrgDb= 'org.Hs.eg.db', 
                                       keyType = 'SYMBOL', 
                                       ont = 'BP',
                                       universe = cts_filt$Genesymbol)

#symbol~TP: 3 GOs found :P
#symbol~TP+reg: No enrichment found in any of gene cluster

#input = genelist_forGO %>% filter(stringr::str_detect(TP, "H9" or "G1"))
#symbol~reg: No enrichment found in any of gene cluster

dotplot(compare.formular_all) + 
  #facet_wrap(~ reg) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10))  +
  scale_y_discrete(labels=function(x) str_wrap(x, width=40))

write.table(compare.formular_all@compareClusterResult, 
            file = here("/Volumes/Extreme SSD/PhD/Drylab/Lancaster/output/GO/2024_01_29_GO_atEachTP.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")

#Second try
#Create a emapplot with human-specifically correlated genes at each time point
#PS: most of the correlated genes in human is human-specific :P
genelist_forGO_2 <- data.frame()

for (i in 1:length(List_H)){
  hsp_forGO <- setdiff(List_H[[i]], List_G[[i]]) %>% 
    as.data.frame() %>% 
    mutate(TP = names(List_H)[i]) %>%
    mutate(state = "h-specific") %>%
    mutate(symbol = ID_symbol %>% 
             filter(Geneid %in% setdiff(List_H[[i]], List_G[[i]])) %>% 
             pull(2))
  overlap <- intersect(List_H[[i]], List_G[[i]]) %>%
    as.data.frame() %>% 
    mutate(TP = names(List_H)[i]) %>%
    mutate(state = "overlaps") %>%
    mutate(symbol = ID_symbol %>% 
             filter(Geneid %in% intersect(List_H[[i]], List_G[[i]])) %>% 
             pull(2))

  genelist_forGO_2 <- rbind(genelist_forGO_2, hsp_forGO, overlap)
}
colnames(genelist_forGO_2)[1] <- "gene"

#compareCluster with the data above
compare.formular_all_2 <- compareCluster(symbol~TP+state, 
                                       data =genelist_forGO_2,
                                       pvalueCutoff = 0.05,
                                       fun = 'enrichGO',
                                       OrgDb= 'org.Hs.eg.db', 
                                       keyType = 'SYMBOL', 
                                       ont = 'BP',
                                       universe = cts_filt$Genesymbol)

#symbol~TP: 1 GO at day 0
#symbol~TP+state: 1 GO H9_d0.h-specific
dotplot(compare.formular_all_2)


#Third try
#Create a emapplot with human-specifically correlated genes throughout the whole time
#PS: most of the correlated genes in human is human-specific :P
genelist_forGO_3 <- data.frame()

List_alltime

hsp_forGO <- setdiff(List_alltime[[1]], List_alltime[[2]]) %>% 
    as.data.frame() %>%
    mutate(state = "h-specific") %>%
    mutate(symbol = ID_symbol %>% 
             filter(Geneid %in% setdiff(List_alltime[[1]], List_alltime[[2]])) %>% 
             pull(2))
overlap <- intersect(List_alltime[[1]], List_alltime[[2]]) %>%
    as.data.frame()  %>%
    mutate(state = "overlaps") %>%
    mutate(symbol = ID_symbol %>% 
             filter(Geneid %in% intersect(List_alltime[[1]], List_alltime[[2]])) %>% 
             pull(2))
  
genelist_forGO_3 <- rbind(genelist_forGO_3, hsp_forGO, overlap)

write.table(hsp_forGO, 
            file = here("/Volumes/Extreme SSD/PhD/Drylab/Lancaster/output/correlation_GTF/2024_01_29_hspecific_correlatedgenes_alltime.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")
write.table(overlap, 
            file = here("/Volumes/Extreme SSD/PhD/Drylab/Lancaster/output/correlation_GTF/2024_01_29_overlap_correlatedgenes_alltime.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")
  
colnames(genelist_forGO_3)[1] <- "gene"

#compareCluster with the data above
compare.formular_all_3 <- compareCluster(symbol~state, 
                                         data =genelist_forGO_3,
                                         pvalueCutoff = 0.05,
                                         fun = 'enrichGO',
                                         OrgDb= 'org.Hs.eg.db', 
                                         keyType = 'SYMBOL', 
                                         ont = 'BP',
                                         universe = cts_filt$Genesymbol)

dotplot(compare.formular_all_3)
write.table(compare.formular_all_3@compareClusterResult, 
            file = here("/Volumes/Extreme SSD/PhD/Drylab/Lancaster/output/GO/2024_01_29_dotplot_GO_BP_correlatedgenes_hspecific_alltime.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")

#Forth try
#Create a emapplot with correlated genes that are overlapping between two species, 
#but showing the opposite direction 

PH_NG <- intersect(List_alltime_pn[["pos.human"]], List_alltime_pn[["neg.gorilla"]]) %>% #72 
  as.data.frame() %>%
  mutate(state = "positive in human, negative in gorilla") %>%
  mutate(symbol = ID_symbol %>% 
           filter(Geneid %in% intersect(List_alltime_pn[["pos.human"]], List_alltime_pn[["neg.gorilla"]])) %>% 
           pull(2))

NH_PG <- intersect(List_alltime_pn[["neg.human"]], List_alltime_pn[["pos.gorilla"]]) %>% #77
  as.data.frame() %>%
  mutate(state = "negative in human, positive in gorilla") %>%
  mutate(symbol = ID_symbol %>% 
           filter(Geneid %in% intersect(List_alltime_pn[["neg.human"]], List_alltime_pn[["pos.gorilla"]])) %>% 
           pull(2))

write.table(PH_NG, 
            file = here("/Volumes/Extreme SSD/PhD/Drylab/Lancaster/output/correlation_GTF/2024_01_29_PH_NG_correlatedgenes_alltime_n72.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")
write.table(NH_PG, 
            file = here("/Volumes/Extreme SSD/PhD/Drylab/Lancaster/output/correlation_GTF/2024_01_29_NH_PG_correlatedgenes_alltime_n77.tsv"),
            quote = F, row.names = F, col.names = T, sep = "\t")

genelist_forGO_4 <- rbind(PH_NG, NH_PG)
colnames(genelist_forGO_4)[1] <- "gene"

#compareCluster with the data above
compare.formular_all_4 <- compareCluster(symbol~state, 
                                         data =genelist_forGO_4,
                                         pvalueCutoff = 0.05,
                                         fun = 'enrichGO',
                                         OrgDb= 'org.Hs.eg.db', 
                                         keyType = 'SYMBOL', 
                                         ont = 'BP',
                                         universe = cts_filt$Genesymbol)

dotplot(compare.formular_all_4)

#############################################################################
#GO - enrichGO separately
#############################################################################
#correlated genes at each timepoint -> no significant GOs

go_H_throughout <- enrichGO(gene=ID_symbol %>% 
                    filter(Geneid %in% List_alltime[[1]]) %>% 
                    pull(2),
                  keyType = "SYMBOL", 
                  OrgDb = "org.Hs.eg.db", 
                  ont = "BP", 
                  universe = cts_filt$Genesymbol,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05)



barplot(go_H_throughout %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_H_throughout)

