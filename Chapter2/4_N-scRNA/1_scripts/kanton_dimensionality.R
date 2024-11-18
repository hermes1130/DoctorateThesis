#############################################################################
#Library loading
#############################################################################
library(Seurat)
library(dplyr)
library(stringr)
library(here)

#############################################################################
#Perform linear dimensional reduction
#Determine the dimensionality of the dataset
#############################################################################
### import the integrated data - takes 1 min on my local
integrated <- readRDS("data/kanton_human_chimp_integrated.rds")

### Perform linear dimensional reduction
integrated <- RunPCA(integrated, verbose = FALSE)
#saveRDS(integrated, file = "integrated_PCA.rds")

### Determine the dimensionality of the dataset

#make elblow function
MakeElbow <- function(obj){
  pct <- obj[["pca"]]@stdev / sum(obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  # Elbow plot to visualize 
  ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
    theme_bw()
}

#Export the image
# 1. Open png file
png("elbow.png", width = 1000, height = 1000)
# 2. Create the plot
MakeElbow(integrated)
# 3. Close the file
dev.off()

#############################################################################
#Cluster the cells
#Run non-linear dimensional reduction (UMAP/tSNE)
#PC = 25
#resolution = 0.1/0.2/0.5/1.0/2.0/5.0 separate scripts on allegro
#############################################################################
integrated <- FindNeighbors(integrated, dims = 1:24)
integrated <- FindClusters(integrated, resolution = 1)

#UMAP
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:25)

png("UMAP_25_0.1_origin1.png", width = 1000, height = 1000)
DimPlot(integrated, reduction = "umap", group.by = "origin")
dev.off()

png("UMAP_25_0.1_origin2.png", width = 2000, height = 1000)
DimPlot(integrated, reduction = "umap", split.by = "origin")
dev.off()

png("UMAP_25_0.1_line.png", width = 4000, height = 1000)
DimPlot(integrated, reduction = "umap", split.by = "line")
dev.off()

png("UMAP_25_0.1.png", width = 1000, height = 1000)
DimPlot(integrated, reduction = "umap")
dev.off()

#tsne
integrated <- RunTSNE(integrated, reduction = "pca", dims = 1:25)

png("TSNE_25_0.1_origin1.png", width = 1000, height = 1000)
DimPlot(integrated, reduction = "tsne", group.by = "origin")
dev.off()

png("TSNE_25_0.1_origin2.png", width = 2000, height = 1000)
DimPlot(integrated, reduction = "tsne", split.by = "origin")
dev.off()

png("TSNE_25_0.1_line.png", width = 4000, height = 1000)
DimPlot(integrated, reduction = "tsne", split.by = "line")
dev.off()

png("TSNE_25_0.1.png", width = 1000, height = 1000)
DimPlot(integrated, reduction = "tsne")
dev.off()

#ZEB2 expression
DefaultAssay(integrated) <- "RNA"

png("UMAP_25_0.1_ZEB2.png", width = 1000, height = 1000)
FeaturePlot(integrated, features = "ZEB2", reduction = "umap")
dev.off()

png("TSNE_25_0.1_ZEB2.png", width = 1000, height = 1000)
FeaturePlot(integrated, features = "ZEB2", reduction = "tsne")
dev.off()

#Some additional dimplots
DimPlot(integrated, reduction = "umap", group.by = "origin", cols = c("#ff80b2ff", "#8dd35fff"))

#############################################################################
#Identify conserved cell type markers
#############################################################################
Markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25)
Markers %>% count(cluster)

write.csv(Markers, file ="markergenes_integ.csv", quote = F, row.names = T, col.names = T)

#############################################################################
#Assigning cell type identity to clusters
#############################################################################
library(rPanglaoDB)
BEC <- getMarkers(include = c('PECAM1', 'VWF'), exclude = c('PDPN', 'ACTA2'))
#############################################################################
#First try setup
#FindNeighbors(dims = 24)
#FindClusters(res = 1)
#FindAllMarkers(only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25)
#############################################################################
integrated <- readRDS("output/RDS/kanton_human_chimp_integrated_tsne_umap_24PCs_0.1res.rds")

#Import the marker genes defined by Sabina - Katja suggested to use a database
kanton.markers <- read.csv(here("../../Volumes/Extreme SSD/PhD/Drylab/Kanton/Kanton_supp/Supplementary_Table_4.csv"), header = T)
rownames(kanton.markers) <- kanton.markers$X
kanton.markers <- kanton.markers[,-1]
kanton.markers %>% count(cluster) #23 cell types

kanton.markers_c <- read.csv(here("../../Volumes/Extreme SSD/PhD/Drylab/Kanton/Kanton_supp/Supplementary_Table_6.csv"), header = T)
rownames(kanton.markers_c) <- kanton.markers_c$X
kanton.markers_c <- kanton.markers_c[,-1]
kanton.markers_c %>% count(cluster) #25 cell types

#my.markers <- read.csv(here("output/markergenes_integ_24.csv"), header = T)
#rownames(my.markers) <- my.markers$X
#my.markers <- my.markers[,-1]

#############################################################################
#PanglaoDB - NOPE
#############################################################################
db.markers <- read.table(here("celltypemarker_db/PanglaoDB_markers_27_Mar_2020.tsv"), sep = '\t',header = T)
#ubiquitousness.index: 0 indicates the gene is not expressed in any cell cluster and 1 indicates that the gene is expressed in all cell clusters.
#Sensitivity: shows how frequently the marker is expressed in cells of the particular cell type
#Specificity: shows how frequently the marker is expressed in cells *NOT* of the particular cell type
#Marker count: indicates how many cell types the marker is used in; i.e. how specific the marker is; doens't exist in the downloaded table

#NOPE: Only 10 cell types with 153 marker genes (Xie et al. 2021)
#db.markers2 <- get(load(file='celltypemarker_db/human_pbmc_marker.rda'))

nrow(db.markers) #6366 rows
db.markers <- 
  db.markers %>% filter(str_detect(species, "Hs")) %>% #6034 rows #remove rows only with mouse gene
  filter(organ == "Brain") %>% #1366 rows #filter rows only with brain
  select(official.gene.symbol, cell.type, ubiquitousness.index, sensitivity_human, specificity_human) #%>% #select only relevant cols
#  filter(sensitivity_human > 0.05) #filter rows only above 0.5 sensitivity (> 0.1 has also 132 marker genes anyway)

#variation w/o sensitivity_human filter
db.markers %>% count(cell.type) #31 cell types out of 31 cell types
nrow(db.markers) #1366

#variation w/ sensitivity_human filter
db.markers %>% count(cell.type) #4 cell types out of 31 cell types
nrow(db.markers) #132

MarkerGenes1 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[1]) %>% slice(1:10) %>% pull(-1)
MarkerGenes2 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[2]) %>% slice(1:10) %>% pull(-1)
MarkerGenes3 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[3]) %>% slice(1:10) %>% pull(-1)
MarkerGenes4 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[4]) %>% slice(1:10) %>% pull(-1)
MarkerGenes5 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[5]) %>% slice(1:10) %>% pull(-1)
MarkerGenes6 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[6]) %>% slice(1:10) %>% pull(-1)

getMarkers(include = c('POU5F1', 'LIN28A', 'TERF1')) #without the not founable gene ids

db.markers %>% 
  filter(official.gene.symbol %in% MarkerGenes1) %>% 
  count(cell.type) %>% 
  arrange(desc(n))
db.markers %>% 
  filter(official.gene.symbol %in% MarkerGenes2) %>% 
  count(cell.type) %>% 
  arrange(desc(n))
db.markers %>% 
  filter(official.gene.symbol %in% MarkerGenes3) %>% 
  count(cell.type) %>% 
  arrange(desc(n))
db.markers %>% 
  filter(official.gene.symbol %in% MarkerGenes4) %>% 
  count(cell.type) %>% 
  arrange(desc(n))
db.markers %>% 
  filter(official.gene.symbol %in% MarkerGenes5) %>% 
  count(cell.type) %>% 
  arrange(desc(n))
db.markers %>% 
  filter(official.gene.symbol %in% MarkerGenes6) %>% 
  count(cell.type) %>% 
  arrange(desc(n))

#Summary of the results in the powerpoint scRNA_Kanton
#Marker genes are mostly expressed in more than one particular cluster.
#Marker gene expression is pretty unspecific, e.g. the difference between pct.1 and pct.2 is pretty high, but pct.2 value is never 0 or near 0 
#Each cluster rather represents an unclear subpopulation, e.g. cluster 0 vs 5 (in my opinion, they should be one cluster)
## Next step: PC = 24, resolution = 0.1 or 0.2

#############################################################################
#Second try setup - FINAL DECISION
#FindNeighbors(dims = 24)
#FindClusters(res = 0.1)
#FindAllMarkers(only.pos = TRUE, min.pct = 0.25,  logfc.threshold = 0.25)
#############################################################################
integrated <- readRDS("output/RDS/kanton_human_chimp_integrated_tsne_umap_24PCs_0.1res.rds")

my.markers <- read.csv(here("output/markergenes_integ_24_0.1.csv"), header = T)
my.markers <- read.csv(here("../../Volumes/Extreme SSD/PhD/Drylab/Kanton/output/markerGenes/markergenes_integ_24_0.1.csv"), header = T)
rownames(my.markers) <- my.markers$X
my.markers <- my.markers[,-1]

#number of marker genes in each cluster
my.markers %>% count(cluster)

MarkerGenes1 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[1]) %>% slice(1:10) %>% pull(-1)
MarkerGenes2 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[2]) %>% slice(1:10) %>% pull(-1)
MarkerGenes3 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[3]) %>% slice(1:10) %>% pull(-1)
MarkerGenes4 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[4]) %>% slice(1:10) %>% pull(-1)
MarkerGenes5 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[5]) %>% slice(1:10) %>% pull(-1)
MarkerGenes6 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[6]) %>% slice(1:10) %>% pull(-1)
MarkerGenes7 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[7]) %>% slice(1:10) %>% pull(-1)
MarkerGenes8 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[8]) %>% slice(1:10) %>% pull(-1)
MarkerGenes9 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[9]) %>% slice(1:10) %>% pull(-1)
MarkerGenes10 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[10]) %>% slice(1:10) %>% pull(-1)
MarkerGenes11 <- my.markers %>% filter(cluster %in% unique(my.markers$cluster)[11]) %>% slice(1:10) %>% pull(-1)

FeaturePlot(integrated, features = MarkerGenes1)
FeaturePlot(integrated, features = MarkerGenes2)
FeaturePlot(integrated, features = MarkerGenes3)
FeaturePlot(integrated, features = MarkerGenes4)
FeaturePlot(integrated, features = MarkerGenes5)
FeaturePlot(integrated, features = MarkerGenes6)
FeaturePlot(integrated, features = MarkerGenes7)
FeaturePlot(integrated, features = MarkerGenes8)
FeaturePlot(integrated, features = MarkerGenes9)
FeaturePlot(integrated, features = MarkerGenes10)
FeaturePlot(integrated, features = MarkerGenes11)
FeaturePlot(integrated, features = MarkerGenes12)

path_out <- here("output/figures/")
for (i in 1:length(unique(my.markers$cluster))) {
  png(paste(path_out, '2023_08_08_FeaturePlot_24PCs_0.1res_Cluster', i-1, '_top10markergenes.png'), width = 1000, height = 781)
  FeaturePlot(integrated, features = my.markers %>% filter(cluster %in% unique(my.markers$cluster)[i]) %>% slice(1:10) %>% pull(-1))
  dev.off()
}
#this didn't work

# Single cell heatmap of feature expression
DefaultAssay(integrated) <- "integrated" 
DoHeatmap(integrated, features = c(MarkerGenes1,
                                   MarkerGenes2,
                                   MarkerGenes3,
                                   MarkerGenes4,
                                   MarkerGenes5,
                                   MarkerGenes6,
                                   MarkerGenes7,
                                   MarkerGenes8,
                                   MarkerGenes9,
                                   MarkerGenes10,
                                   MarkerGenes11))

#############################################################################
#Sabina's paper 
#############################################################################
#human + chimp
kanton.markers$species <- "human"
kanton.markers_c$species <- "chimp"
KantonMarker <- rbind(kanton.markers, kanton.markers_c)

#human only
CTmarkers <- read.csv(here("../../Volumes/Extreme SSD/PhD/Drylab/Kanton/output/markerGenes/markers_fromPaper.txt"), header = T, sep = "\t")

#top 10 or 100 markers with lowest p-value and this will be compared with the imported table
top10 <- my.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
top100 <- my.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC)


#a custom function to search the most matching cell type to our markers. 
#For each cluster, the function will check how often the top 10 markers match to certain cell types 
#based on the imported table. 
#The mostly matched cell type for each cluster will be returned.
CTretrieve <- function(x){
  cluster <- c()
  for(i in 1:length(unique(x$cluster))){
    CTtmp <- x %>% filter(cluster %in% unique(x$cluster)[i]) %>% pull(-1)
    cluster <- c(cluster, CTMarker.c %>% 
                   filter(gene %in% CTtmp) %>% 
                   count(cluster) %>% 
                   arrange(desc(n)) %>% 
                   slice(1) %>% 
                   pull(1))
  }
  return(cluster)
}

#only human no merging
#function with CTmarkers
CTretrieve(top10) #overlapping cell types
CTretrieve(top100) # no overlapping cell types

#only chimp no merging
#function with kanton.markers_c
CTretrieve(top10) #overlapping cell types
CTretrieve(top100) # no overlapping cell types

#human+chimp no merging
#function with KantonMarker
CTretrieve(top10) #overlapping cell types
CTretrieve(top100) #overlapping cell types

#merge the cell type 1 and 2 (or even more) together
CTmarkers %>% count(cluster) #23 cell types
mergCTmarkers <- CTmarkers
mergCTmarkers$cluster <- gsub('[[:digit:]]+', '', CTmarkers$cluster)
mergCTmarkers %>% count(cluster) #15 cell types

#only human with merging
#function with mergCTmarkers
CTretrieve(top10) #overlapping cell types
CTretrieve(top100) #more overlapping cell types

kanton.markers_c %>% count(cluster) #25 cell types
CTMarker.c <- kanton.markers_c
CTMarker.c$cluster <- gsub('[[:digit:]]+', '', kanton.markers_c$cluster)
CTMarker.c %>% count(cluster) #14 cell types

#only chimp with merging
#function with CTMarker.c
CTretrieve(top10) #overlapping cell types
CTretrieve(top100) #overlapping cell types

KantonMarker %>% count(cluster) #41 cell types
mergKantonMarker <- KantonMarker
mergKantonMarker$cluster <- gsub('[[:digit:]]+', '', KantonMarker$cluster)
mergKantonMarker %>% count(cluster) #25 cell types

#human+chimp with merging
#function with mergKantonMarker
CTretrieve(top10) #overlapping cell types
CTretrieve(top100) #more overlapping cell types

#After discussion with Katja and Vladi, we decide to stick with only many tries with Sabina's cell types
#First check if cluster 1 & 7 and cluster 
c(GetGeneName.pos.ct(cor_H, 1), GetGeneName.neg.ct(cor_H, 1))





  
  
  

