library(Seurat)
library(dplyr)
library(stringr)
library(tidyverse)


#############################################################################
#Assigning cell type identity to clusters - not now
#############################################################################
#list of cell types containing ZEB2 expression
#CT_human <- c("cortical neurons 1",  "ventral progenitors and neurons 3","NSC", "G2M/S ventral progenitors and IPs", 
#              "NSC/radial glia","gliogenic/outer radial glia", "IPs and early cortical neurons", "ventral progenitors and neurons 2",
#              "G2M/S NPCs", "cortical neurons 2", "cortical neurons 3",
#              "ventral progenitors and neurons 4", "mesenchymal-like cells", "G2M/S dorsal progenitors 3")
#CT_chimp <- c("cortical neurons 2","upper layer neurons 2","upper layer neurons 3",
#              "deep layer neurons 2", "cortical neurons 1","deep layer neurons 1",
#              "ventral forebrain progenitors and neurons", "G2M/S dorsal progenitors","mesenchymal-like cells","unknown")

#############################################################################
#running the correlation test on the server - including all cells
#############################################################################
# Subset metadata
meta <- as.data.frame(integrated@meta.data %>% select(seurat_clusters))
colnames(meta) <- "CT"
dim(meta)

ZEB2_row <- which(rownames(integrated@assays$SCT@data[,1:nrow(meta)]) =="ZEB2")
ZEB2_row

#Correlation test within human samples
tmp_cor <- data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=1))
rownames(tmp_cor) <- rownames(integrated@assays$SCT@data)
colnames(tmp_cor) <- "cor"
dim(tmp_cor)

cor_H <-data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=length(table(Idents(integrated)))))
rownames(cor_H) <- rownames(integrated@assays$SCT@data)
colnames(cor_H) <- names(table(Idents(integrated)))
dim(cor_H)
colnames(cor_H)

for (i in 1:length(names(table(Idents(integrated))))){
  tmp_stage <- as.matrix(integrated@assays$SCT@data[, rownames(meta %>% dplyr::filter(stringr::str_detect(rownames(meta), "Human")) %>% dplyr::filter(stringr::str_detect(CT, names(table(Idents(integrated)))[i])))])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- cor.test(as.numeric(tmp_stage[ZEB2_row,]), as.numeric(tmp_stage[j,]),method ="pearson")$p.value
    tmp_cor[is.na(tmp_cor)] <- 1
    tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
    if (tmp_cor[j,] < 0.05){
      cor_H[j, names(table(Idents(integrated)))[i]] <- cor.test(as.numeric(tmp_stage[ZEB2_row,]), as.numeric(tmp_stage[j,]),method ="pearson")$estimate
    } 
  }
}

saveRDS(cor_H, "cor_H.rds")

#Correlation test within chimpanzee samples
tmp_cor <- data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=1))
rownames(tmp_cor) <- rownames(integrated@assays$SCT@data)
colnames(tmp_cor) <- "cor"
dim(tmp_cor)

cor_C <-data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=length(table(Idents(integrated)))))
rownames(cor_C) <- rownames(integrated@assays$SCT@data)
colnames(cor_C) <- names(table(Idents(integrated)))
dim(cor_C)
colnames(cor_C)

for (i in 1:length(names(table(Idents(integrated))))){
  tmp_stage <- as.matrix(integrated@assays$SCT@data[, rownames(meta %>% dplyr::filter(stringr::str_detect(rownames(meta), "Chimp")) %>% dplyr::filter(stringr::str_detect(CT, names(table(Idents(integrated)))[i])))])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- cor.test(as.numeric(tmp_stage[ZEB2_row,]), as.numeric(tmp_stage[j,]),method ="pearson")$p.value
    tmp_cor[is.na(tmp_cor)] <- 1
    tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
    if (tmp_cor[j,] < 0.05){
      cor_C[j, names(table(Idents(integrated)))[i]] <- cor.test(as.numeric(tmp_stage[ZEB2_row,]), as.numeric(tmp_stage[j,]),method ="pearson")$estimate
    } 
  }
}

saveRDS(cor_C, "cor_C.rds")

#############################################################################
#How many genes are expressing ZEB2 in each cluster*species 
#############################################################################
#Run this part on the server
# Subset metadata
meta <- as.data.frame(integrated@meta.data %>% select(seurat_clusters))
colnames(meta) <- "CT"
dim(meta)

ZEB2_row <- which(rownames(integrated@assays$SCT@data[,1:nrow(meta)]) =="ZEB2")
ZEB2_row

#extract the number of ZEB2 expressing gene in human samples (> 0.25)
n <- c()
for (i in 1:length(names(table(Idents(integrated))))){
  tmp_stage <-  as.matrix(integrated@assays$SCT@data[, rownames(meta %>% dplyr::filter(stringr::str_detect(rownames(meta), "Human")) %>% dplyr::filter(stringr::str_detect(CT, names(table(Idents(integrated)))[i])))])
  n <- c(n, sum(tmp_stage[ZEB2_row,] > 0.25))
}

saveRDS(n, "ZEB2expressingGenes_eachCluster_0.25_Human.rds")

#extract the number of ZEB2 expressing gene in chimp samples (> 0.25)
m <- c()
for (i in 1:length(names(table(Idents(integrated))))){
  tmp_stage <-  as.matrix(integrated@assays$SCT@data[, rownames(meta %>% dplyr::filter(stringr::str_detect(rownames(meta), "Chimp")) %>% dplyr::filter(stringr::str_detect(CT, names(table(Idents(integrated)))[i])))])
  m <- c(m, sum(tmp_stage[ZEB2_row,] > 0.25))
}

saveRDS(m, "ZEB2expressingGenes_eachCluster_0.25_Chimp.rds")

rownames(integrated@assays$SCT@data[rowSums(integrated@assays$SCT@data)>0,]) #29338
nrow(integrated@assays$SCT@data) #29339

AllGenes <- rownames(integrated@assays$SCT@data) #BackgroundGenes_forGO_allexpressedgenes_afterFilter.tsv
write.table(AllGenes, file = paste0(path_out, "BackgroundGenes_forGO_allgenes.tsv"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

#############################################################################
#running the correlation test on the server - including only ZEB2 expressing cells
#############################################################################
#ZEB2 expressing cells are defined as non-zero normalized count
# Subset metadata
meta <- as.data.frame(integrated@meta.data %>% select(seurat_clusters))
colnames(meta) <- "CT"
dim(meta)

ZEB2_row <- which(rownames(integrated@assays$SCT@data[,1:nrow(meta)]) =="ZEB2")
ZEB2_row

#Correlation test within human samples
tmp_cor <- data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=1))
rownames(tmp_cor) <- rownames(integrated@assays$SCT@data)
colnames(tmp_cor) <- "cor"
dim(tmp_cor)

cor_H <-data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=length(table(Idents(integrated)))))
rownames(cor_H) <- rownames(integrated@assays$SCT@data)
colnames(cor_H) <- names(table(Idents(integrated)))
dim(cor_H)
colnames(cor_H)

for (i in 1:length(names(table(Idents(integrated))))){
  tmp_stage <- as.matrix(integrated@assays$SCT@data[, rownames(meta %>% 
                                                                 dplyr::filter(stringr::str_detect(rownames(meta), "Human")) %>% 
                                                                 dplyr::filter(stringr::str_detect(CT, names(table(Idents(integrated)))[i])))])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- cor.test(as.numeric(tmp_stage[ZEB2_row, which(tmp_stage[ZEB2_row,] > 0)]), 
                            as.numeric(tmp_stage[j,which(tmp_stage[ZEB2_row,] > 0)]),
                            method ="pearson")$p.value
    tmp_cor[is.na(tmp_cor)] <- 1
    tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
    if (tmp_cor[j,] < 0.05){
      cor_H[j, names(table(Idents(integrated)))[i]] <- cor.test(as.numeric(tmp_stage[ZEB2_row,which(tmp_stage[ZEB2_row,] > 0)]), 
                                                                as.numeric(tmp_stage[j,which(tmp_stage[ZEB2_row,] > 0)]),
                                                                method ="pearson")$estimate
    } 
  }
}

saveRDS(cor_H, "cor_H_onlyZEB2expressingcells.rds")

#Correlation test within chimp samples
tmp_cor <- data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=1))
rownames(tmp_cor) <- rownames(integrated@assays$SCT@data)
colnames(tmp_cor) <- "cor"
dim(tmp_cor)

cor_C <-data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=length(table(Idents(integrated)))))
rownames(cor_C) <- rownames(integrated@assays$SCT@data)
colnames(cor_C) <- names(table(Idents(integrated)))
dim(cor_C)
colnames(cor_C)

#slight change of the i definition because cluster 10 has only 1 ZEB2 expressing cell
#smaller than 3 samples - no correlation test possible
for (i in 2:length(names(table(Idents(integrated))))-1){
  tmp_stage <- as.matrix(integrated@assays$SCT@data[, rownames(meta %>% 
                                                                 dplyr::filter(stringr::str_detect(rownames(meta), "Chimp")) %>% 
                                                                 dplyr::filter(stringr::str_detect(CT, names(table(Idents(integrated)))[i])))])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- cor.test(as.numeric(tmp_stage[ZEB2_row, which(tmp_stage[ZEB2_row,] > 0)]), 
                            as.numeric(tmp_stage[j,which(tmp_stage[ZEB2_row,] > 0)]),
                            method ="pearson")$p.value
    tmp_cor[is.na(tmp_cor)] <- 1
    tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
    if (tmp_cor[j,] < 0.05){
      cor_C[j, names(table(Idents(integrated)))[i]] <- cor.test(as.numeric(tmp_stage[ZEB2_row,which(tmp_stage[ZEB2_row,] > 0)]), 
                                                                as.numeric(tmp_stage[j,which(tmp_stage[ZEB2_row,] > 0)]),
                                                                method ="pearson")$estimate
    } 
  }
}

saveRDS(cor_C, "cor_C_onlyZEB2expressingcells.rds")

#############################################################################
#running the correlation test on the server - including only ZEB2 expressing cells
#merged and assigned CTs
#############################################################################
#ZEB2 expressing cells are defined as non-zero normalized count
# Subset metadata
meta <- as.data.frame(integrated@meta.data %>% select(seurat_clusters))
colnames(meta) <- "CT"
dim(meta)
str(meta)
#convert CT variable to character
meta$CT <- as.character(meta$CT)
str(meta)

#List of CTs according to the cluster number
CTList_merg <- c("cortical neurons",
                 "GM/S dorsal progenitors/radial glia",
                 "stem cells",
                 "ventral progenitors and neurons",
                 "neuroectodermal-like cells",
                 "midbrain/hindbrain",
                 "NSC/radial glia",
                 "GM/S dorsal progenitors/radial glia",
                 "NSCs",
                 "mesenchymal-like cells",
                 "stem cells")

for(i in 1:length(names(table(Idents(integrated))))){
  meta[meta == names(table(Idents(integrated)))[i]] <- CTList_merg[i]
}

ZEB2_row <- which(rownames(integrated@assays$SCT@data[,1:nrow(meta)]) =="ZEB2")
ZEB2_row

#UMAP Dimplot with merged CT
integrated@meta.data$assignedCT <- meta$CT
DimPlot(integrated, reduction = "umap", group.by = "assignedCT")
#UMAP Dimplot with "stage" 
unique(integrated@meta.data$stage_group)
DimPlot(integrated, reduction = "umap", group.by = "stage_group")
#UMAP Dimplot split by species
DimPlot(integrated, reduction = "umap", split.by = "species")

#UMAP FeaturePlot with ZEB2 expressing cells
FeaturePlot(integrated, reduction = "umap", features = "ZEB2")
#UMAP FeaturePlot with ZEB2 expressing cells + split by species
FeaturePlot(integrated, reduction = "umap", features = "ZEB2", split.by = "species")


#Correlation test within human samples
tmp_cor <- data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=1))
rownames(tmp_cor) <- rownames(integrated@assays$SCT@data)
colnames(tmp_cor) <- "cor"
dim(tmp_cor)

cor_H <-data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=length(unique(meta$CT))))
rownames(cor_H) <- rownames(integrated@assays$SCT@data)
colnames(cor_H) <- unique(meta$CT)
dim(cor_H)
colnames(cor_H)

for (i in 1:length(unique(meta$CT))){
  tmp_stage <- as.matrix(integrated@assays$SCT@data[, rownames(meta %>% 
                                                                 dplyr::filter(stringr::str_detect(rownames(meta), "Human")) %>% 
                                                                 dplyr::filter(stringr::str_detect(CT, unique(meta$CT)[i])))])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- cor.test(as.numeric(tmp_stage[ZEB2_row, which(tmp_stage[ZEB2_row,] > 0)]), 
                            as.numeric(tmp_stage[j,which(tmp_stage[ZEB2_row,] > 0)]),
                            method ="pearson")$p.value
    tmp_cor[is.na(tmp_cor)] <- 1
    tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
    if (tmp_cor[j,] < 0.05){
      cor_H[j, unique(meta$CT)[i]] <- cor.test(as.numeric(tmp_stage[ZEB2_row,which(tmp_stage[ZEB2_row,] > 0)]), 
                                                                as.numeric(tmp_stage[j,which(tmp_stage[ZEB2_row,] > 0)]),
                                                                method ="pearson")$estimate
    } 
  }
}

saveRDS(cor_H, "cor_H_onlyZEB2expressingcells_mergednassignedCT.rds")

#Correlation test within chimp samples
tmp_cor <- data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=1))
rownames(tmp_cor) <- rownames(integrated@assays$SCT@data)
colnames(tmp_cor) <- "cor"
dim(tmp_cor)

cor_C <-data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=length(unique(meta$CT))))
rownames(cor_C) <- rownames(integrated@assays$SCT@data)
colnames(cor_C) <- unique(meta$CT)
dim(cor_C)
colnames(cor_C)

for (i in 1:length(unique(meta$CT))){
  tmp_stage <- as.matrix(integrated@assays$SCT@data[, rownames(meta %>% 
                                                                 dplyr::filter(stringr::str_detect(rownames(meta), "Chimp")) %>% 
                                                                 dplyr::filter(stringr::str_detect(CT, unique(meta$CT)[i])))])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- cor.test(as.numeric(tmp_stage[ZEB2_row, which(tmp_stage[ZEB2_row,] > 0)]), 
                            as.numeric(tmp_stage[j,which(tmp_stage[ZEB2_row,] > 0)]),
                            method ="pearson")$p.value
    tmp_cor[is.na(tmp_cor)] <- 1
    tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
    if (tmp_cor[j,] < 0.05){
      cor_C[j, unique(meta$CT)[i]] <- cor.test(as.numeric(tmp_stage[ZEB2_row,which(tmp_stage[ZEB2_row,] > 0)]), 
                                                                as.numeric(tmp_stage[j,which(tmp_stage[ZEB2_row,] > 0)]),
                                                                method ="pearson")$estimate
    } 
  }
}

saveRDS(cor_C, "cor_C_onlyZEB2expressingcells_mergednassignedCT.rds")

#############################################################################
#repeating the test
#############################################################################
#rationale of this strategy is that the number of the ZEB2 expressing cells differs
#very much between two species, and now I want to test multiple times with the same
#number of the cells 
#test only "cortical neurons", "GM/S dorsal progenitors/radial glia",
#"mesenchymal-like cells", "ventral progenitors and neurons"


#ZEB2 expressing cells are defined as non-zero normalized count
# Subset metadata
meta <- as.data.frame(integrated@meta.data %>% select(seurat_clusters))
colnames(meta) <- "CT"
dim(meta)
str(meta)
#convert CT variable to character
meta$CT <- as.character(meta$CT)
str(meta)

#List of CTs according to the cluster number
CTList_merg <- c("cortical neurons",
                 "GM/S dorsal progenitors/radial glia",
                 "stem cells",
                 "ventral progenitors and neurons",
                 "neuroectodermal-like cells",
                 "midbrain/hindbrain",
                 "NSC/radial glia",
                 "GM/S dorsal progenitors/radial glia",
                 "NSCs",
                 "mesenchymal-like cells",
                 "stem cells")

for(i in 1:length(names(table(Idents(integrated))))){
  meta[meta == names(table(Idents(integrated)))[i]] <- CTList_merg[i]
}
unique(meta$CT)

ZEB2_row <- which(rownames(integrated@assays$SCT@data[,1:nrow(meta)]) =="ZEB2")
ZEB2_row

tmp_cor <- data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol=1))
rownames(tmp_cor) <- rownames(integrated@assays$SCT@data)
colnames(tmp_cor) <- "cor"
dim(tmp_cor)

for (i in 1:length(unique(meta$CT))){
  tmp_stage_h <- as.matrix(integrated@assays$SCT@data[, rownames(meta %>% 
                                                                 dplyr::filter(stringr::str_detect(rownames(meta), "Human")) %>% 
                                                                 dplyr::filter(stringr::str_detect(CT, unique(meta$CT)[i])))])
  tmp_stage_c <- as.matrix(integrated@assays$SCT@data[, rownames(meta %>% 
                                                                   dplyr::filter(stringr::str_detect(rownames(meta), "Chimp")) %>% 
                                                                   dplyr::filter(stringr::str_detect(CT, unique(meta$CT)[i])))])
  n_ZEB2_h <- length(tmp_stage_h[1, which(tmp_stage_h[1,] > 0)])
  n_ZEB2_c <- length(tmp_stage_c[1, which(tmp_stage_c[1,] > 0)])
  
  if (n_ZEB2_h > n_ZEB2_c){
    #number of zeb2 expressing cells is higher in human
    for (j in 1:trunc(n_ZEB2_h/n_ZEB2_c)-1){
      #zeb2 expressing cells only in human
      tmp_stage_h_zeb2 <- tmp_stage_h[, which(tmp_stage_h[ZEB2_row,] > 0 )]
      #zeb2 expressing cells in human with the same number of chimp sample
      tmp_stage_h_zeb2_div <- tmp_stage_h_zeb2[,(n_ZEB2_h*j)+1:n_ZEB2_h*(j+1)]
      ncol(tmp_stage_h_zeb2_div)
      tmp_stage_h_zeb2_tail <- tmp_stage_h_zeb2[,c(ncol(tmp_stage_h_zeb2)-n_ZEB2_h+1):ncol(tmp_stage_h_zeb2)]
      ncol(tmp_stage_c_zeb2_tail)
      
      #data prep
      
      cor_h_div <-data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol= trunc(n_ZEB2_h/n_ZEB2_c)+1))
      rownames(cor_h_div) <- rownames(integrated@assays$SCT@data)
      colnames(cor_h_div) <- c(1:trunc(n_ZEB2_h/n_ZEB2_c)+1)
      dim(cor_h_div)
      colnames(cor_h_div)
      
      #correlation test
      
      for (l in 1:nrow(tmp_stage_h_zeb2_div)){
        tmp_cor[l,] <- cor.test(as.numeric(tmp_stage_h_zeb2_div[ZEB2_row, ]), 
                                as.numeric(tmp_stage_h_zeb2_div[l,]),
                                method ="pearson")$p.value
        tmp_cor[is.na(tmp_cor)] <- 1
        tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
        if (tmp_cor[l,] < 0.05){
          cor_h_div[l, j+1] <- cor.test(as.numeric(tmp_stage_h_zeb2_div[ZEB2_row,]), 
                                        as.numeric(tmp_stage_h_zeb2_div[l,]),
                                        method ="pearson")$estimate
        } 
      }
      
      for (m in 1:nrow(tmp_stage_h_zeb2_tail)){
        tmp_cor[m,] <- cor.test(as.numeric(tmp_stage_h_zeb2_tail[ZEB2_row, ]), 
                                as.numeric(tmp_stage_h_zeb2_tail[m,]),
                                method ="pearson")$p.value
        tmp_cor[is.na(tmp_cor)] <- 1
        tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
        if (tmp_cor[m,] < 0.05){
          cor_h_div[m, trunc(n_ZEB2_c/n_ZEB2_h)+1] <- cor.test(as.numeric(tmp_stage_h_zeb2_tail[ZEB2_row,]), 
                                                               as.numeric(tmp_stage_h_zeb2_tail[m,]),
                                                               method ="pearson")$estimate
        } 
      }
      
      saveRDS(cor_h_div, paste('cor_H_only_onlyZEB2expressingcells_mergednassignedCT_div_',
                               i,
                               '.rds'))  
      
    }
  } else {
    #number of zeb2 expressing cells is higher in chimp
    for (j in 0:trunc(n_ZEB2_c/n_ZEB2_h)-1){
      #zeb2 expressing cells only in chimp
      tmp_stage_c_zeb2 <- tmp_stage_c[, which(tmp_stage_c[ZEB2_row,] > 0 )]
      #zeb2 expressing cells in chimp with the same number of human sample
      tmp_stage_c_zeb2_div <- tmp_stage_c_zeb2[,(n_ZEB2_h*j)+1:n_ZEB2_h*(j+1)]
      ncol(tmp_stage_c_zeb2_div)
      tmp_stage_c_zeb2_tail <- tmp_stage_c_zeb2[,c(ncol(tmp_stage_c_zeb2)-n_ZEB2_h+1):ncol(tmp_stage_c_zeb2)]
      ncol(tmp_stage_c_zeb2_tail)
      
      #data prep
      
      cor_c_div <-data.frame(matrix(, nrow=nrow(integrated@assays$SCT@data), ncol= trunc(n_ZEB2_c/n_ZEB2_h)+1))
      rownames(cor_c_div) <- rownames(integrated@assays$SCT@data)
      colnames(cor_c_div) <- c(1:trunc(n_ZEB2_c/n_ZEB2_h)+1)
      dim(cor_c_div)
      colnames(cor_c_div)
      
      #correlation test
      
      for (l in 1:nrow(tmp_stage_c_zeb2_div)){
        tmp_cor[l,] <- cor.test(as.numeric(tmp_stage_c_zeb2_div[ZEB2_row, ]), 
                                as.numeric(tmp_stage_c_zeb2_div[l,]),
                                method ="pearson")$p.value
        tmp_cor[is.na(tmp_cor)] <- 1
        tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
        if (tmp_cor[l,] < 0.05){
          cor_c_div[l, j+1] <- cor.test(as.numeric(tmp_stage_c_zeb2_div[ZEB2_row,]), 
                                        as.numeric(tmp_stage_c_zeb2_div[l,]),
                                        method ="pearson")$estimate
        } 
      }
      
      for (m in 1:nrow(tmp_stage_c_zeb2_tail)){
        tmp_cor[m,] <- cor.test(as.numeric(tmp_stage_c_zeb2_tail[ZEB2_row, ]), 
                                as.numeric(tmp_stage_c_zeb2_tail[m,]),
                                method ="pearson")$p.value
        tmp_cor[is.na(tmp_cor)] <- 1
        tmp_cor[,1] <- p.adjust(tmp_cor[,1], method = "BH")
        if (tmp_cor[m,] < 0.05){
          cor_c_div[m, trunc(n_ZEB2_c/n_ZEB2_h)+1] <- cor.test(as.numeric(tmp_stage_c_zeb2_tail[ZEB2_row,]), 
                                                               as.numeric(tmp_stage_c_zeb2_tail[m,]),
                                                               method ="pearson")$estimate
        } 
      }
      
      saveRDS(cor_c_div, paste('cor_C_only_onlyZEB2expressingcells_mergednassignedCT_div_',
                               i,
                               '.rds'))     
      
    }
  }
  
}




