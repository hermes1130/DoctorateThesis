#############################################################################
#Library loading
#############################################################################
library(Seurat)
library(dplyr)
#############################################################################
#Data upload, containing cell ID and cell type information
#############################################################################
tc_chimp_sc <- readRDS(here("data/timecourse_chimp_singlecells_consensusGenome.rds"))
tc_human_sc <- readRDS(here("data/timecourse_human_singlecells_GRCh38.rds"))
#############################################################################
#how many cell lines?
table(tc_human_sc$line) #h409B2:20272     H9:23226 
table(tc_chimp_sc$line) #Sandra:31611     SandraA:5273 

# split the dataset into a list of two seurat objects by cell lines
Split_tc_human <- SplitObject(tc_human_sc, split.by = "line")
Split_tc_chimp <- SplitObject(tc_chimp_sc, split.by = "line")

tc_human_sc_sct <- SCTransform(tc_human_sc, vars.to.regress = "percent.mito", verbose = FALSE)

#############################################################################
#Server script
#############################################################################
#Data import
tc_human <- readRDS(file = "/data/scratch/hermes1130/timecourse_human_singlecells_GRCh38.rds")
tc_chimp <- readRDS(file = "/data/scratch/hermes1130/timecourse_chimp_singlecells_consensusGenome.rds")

#SCtransform
tc_human_sct <- SCTransform(tc_human, vars.to.regress = "percent.mito", verbose = FALSE)
tc_chimp_sct <- SCTransform(tc_chimp, vars.to.regress = "percent.mito", verbose = FALSE)

#Save RDS - SCtransform
saveRDS(tc_human_sct, file = "tc_human_sct.rds")
saveRDS(tc_chimp_sct, file = "tc_chimp_sct.rds")

#Merging Two Seurat Objects
merged <- merge(tc_human_sct, y = tc_chimp_sct, add.cell.ids = c("Human", "Chimp"), project = "HSvsPT", merge.data = TRUE)
merged
head(colnames(merged))

#define variable features
VariableFeatures(merged[["SCT"]]) <- rownames(merged[["SCT"]]@scale.data)
head(VariableFeatures(merged[["SCT"]]))

#add a new column
merged$origin <- str_detect(rownames(merged[[]]), "Human")
#TRUE == Human
#FALSE == Chimp
merged <- readRDS(file = "data/merged.rds")
head(merged@meta.data, 5)

#split the dataset by "origin"
splitted <- SplitObject(merged, split.by = "origin")
splitted

#select integration features
features <- SelectIntegrationFeatures(object.list = splitted, nfeatures = 3000)
head(features)

#prep for SCT integration
splitted_prep <- PrepSCTIntegration(object.list = merged, anchor.features = features)
splitted_prep

#find anchors
anchors <- FindIntegrationAnchors(object.list = splitted_prep, normalization.method = "SCT", anchor.features = features)
anchors

#perform integration
integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
integrated

#Save RDS
saveRDS(integrated, file = "kanton_human_chimp_integrated.rds")

