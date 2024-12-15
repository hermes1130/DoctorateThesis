#cor_H <- readRDS(here("output/RDS/cor_H_final_HC_integSCT.rds"))
#cor_C <- readRDS(here("output/RDS/cor_C_final_HC_integSCT.rds"))
cor_H <- readRDS(file = "../../Volumes/Extreme SSD/PhD/Drylab/Kanton/output/RDS/cor_H_onlyZEB2expressingcells.rds")
cor_C <- readRDS(here("output/RDS/cor_C_onlyZEB2expressingcells.rds"))
dim(cor_H)
dim(cor_C)

#remove ZEB2 which is counted as one of the correlated genes
cor_H <- cor_H[-which(rownames(cor_H) == "ZEB2"),]
cor_C <- cor_C[-which(rownames(cor_C) == "ZEB2"),]
dim(cor_H)
dim(cor_C)

#############################################################################
#Plotting the number of correlated genes with ZEB2 in each CT
#############################################################################

nr_corr <- data.frame(matrix(, nrow=c(ncol(cor_H)+ncol(cor_C))*2, ncol=4))
colnames(nr_corr) <- c("CT", "correlation", "species","n")
nr_corr[,1] <- c(rep(colnames(cor_H), 2), rep(colnames(cor_C),2))
nr_corr[,2] <- c(rep("pos", ncol(cor_H)), rep("neg",ncol(cor_H)), rep("pos", ncol(cor_C)), rep("neg", ncol(cor_C)))
nr_corr[,3] <- c(rep("Human", ncol(cor_H)*2), rep("Chimp", ncol(cor_C)*2))

for (i in 1:nrow(nr_corr)){
  if (i < ncol(cor_H)*2+1) {
    if (i < ncol(cor_H)+1){
      nr_corr[i,4] <- as.numeric(length(rownames(cor_H)[which(sign(cor_H[,i]) == 1)]))
    } else {
      nr_corr[i,4] <- as.numeric(length(rownames(cor_H)[which(sign(cor_H[,i-11])==-1)]))
    }
  } else {
    if (i < ncol(cor_H)*2+ncol(cor_C)+1) {
      nr_corr[i,4] <- as.numeric(length(rownames(cor_C)[which(sign(cor_C[,i-22])==1)]))
    } else {
      nr_corr[i,4] <- as.numeric(length(rownames(cor_C)[which(sign(cor_C[,i-33])==-1)]))
    }
  }
}


#change the nr of corr into numeric
nr_corr$n <- as.numeric(nr_corr$n)

path_out <- here("output/correlation/")
write.table(nr_corr, 
            file = "output/correlation/overview_tbl_nr_correlation_eachcluster_INTEG_onlyZEB2expressingcells.tsv", 
            row.names = F, col.names = T, quote = F,sep = "\t")


library(plyr)
# Sort by timpoint and species
nr_corr_sort <- arrange(nr_corr, CT, species)
# Calculate the cumulative sum of n for each correlation
nr_corr_sum <- ddply(nr_corr_sort, c("CT", "species"),
                     transform, label_ypos=cumsum(n))

#calculate the log values of counts
nr_corr_sum$log_n <- log(nr_corr_sum$n)
#change the infinite values to zero
nr_corr_sum$log_n[which(!is.finite(nr_corr_sum$log_n))] <- 0
#calculate the log values of cumulative sums
nr_corr_sum$log_label_ypos <- log(nr_corr_sum$label_ypos)

#order the levels
nr_corr_sum$CT <- factor(nr_corr_sum$CT, levels =colnames(cor_H))
nr_corr_sum$species <- factor(nr_corr_sum$species, levels =c("Human", "Chimp"))

#barplot with number of correlated genes
ggplot(data=filter(nr_corr_sum), aes(x=CT, y=n, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.3, color="white", size=2.5)+
  scale_y_continuous(limits=c(0, 1300), breaks = c(0, 100, 250, 500, 1000))+
  facet_wrap( ~ species) +
  theme(axis.title = element_blank()) +
  labs(title="Number of correlated genes with ZEB2 expression", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2") + 
  scale_y_break(c(260, 400, 500, 800), scales = 0.2) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))

#percent stacked barplot showing number of correlated genes per cluser
my_cols <- c('0'='#F8766D','1'='#DB8E00','2'='#AEA200','3'='#64B200','4'='#00BD5C',
             '5'='#00C1A7','6'='#00BADE','7'='#00A6FF','8'='#B385FF','9'='#EF67EB',
             '10'='#FF63B6')

ggplot(data=filter(nr_corr_sum), aes(x=species, y=n, fill=CT)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.title = element_blank()) +
  labs(title="Proportion of correlated genes in each cluster", x="", y="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))

ggplot(data=filter(nr_corr_sum), aes(x=correlation, y=n, fill=CT)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.title = element_blank()) +
  facet_wrap( ~ species) +
  labs(title="Proportion of correlated genes in each cluster", x="", y="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))

#############################################################################
ct_h <- ggplot(data=filter(filter(nr_corr_sum), species == "human"), aes(x=CT, y=n, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.3, color="white", size=2.5)+
  scale_y_continuous(limits=c(0, 9000), breaks=c(0, 4000, 5000, 6000, 7000, 8000))+
  scale_y_break(c(10, 15), scales = 2) +
  scale_y_break(c(30, 50), scales = 3) +
  scale_y_break(c(250, 800), scales = 1.5) +
  facet_wrap( ~ species) +
  theme(axis.title = element_blank()) +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))

ct_c <-ggplot(data=filter(filter(nr_corr_sum, !grepl("cortical neurons 1|cortical neurons 2|mesenchymal-like cells", CT)), species == "chimpanzee"), aes(x=CT, y=n, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.3, color="white", size=2.5)+
  scale_y_continuous(limits=c(0, 690), breaks=c(0, 5, 50, 100, 200, 300, 600))+
  scale_y_break(c(3, 45), scales = 5) +
  scale_y_break(c(290, 450), scales = 1.5) +
  facet_wrap( ~ species) +
  theme(axis.title = element_blank()) +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))


ct3 <- ggplot(data=filter(nr_corr_sum, grepl("cortical neurons 1|cortical neurons 2|mesenchymal-like cells", CT)), 
       aes(x=species, y=n, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.5, color="white", size=2.5)+
  scale_y_continuous(limits=c(0, 250), breaks=c(0, 5, 25, 50, 100, 250))+
  scale_y_break(c(6, 16), scales = 1.5) +
  scale_y_break(c(55, 90), scales = 1.5) +
  facet_wrap( ~ CT) +
  theme(axis.title = element_blank()) +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))

#labs(title="Number of correlated genes with ZEB2 expression", x="") 
ggarrange(ct3, ct_h, ct_c, nrow = 1)

#Size of ggplot above: 950x530

#############################################################################
#ZEB2 expressing samples per cluster*species
#############################################################################
#ZEB2 expression per cluster
VlnPlot(integrated, features = "ZEB2", assay = "RNA") + geom_hline(yintercept = 0.5, color="red", linetype="dashed")
VlnPlot(integrated, features = "ZEB2", assay = "SCT") + geom_hline(yintercept = 0.25, color="red", linetype="dashed")
#Density plot
integrated[[]] %>%
  ggplot(aes(x = nCount_SCT + 1)) + 
  geom_density(color = "gray80", linetype = 2, size = 1.5)+ 
  geom_density(aes(color = seurat_clusters)) +
  scale_x_log10() +
  theme_bw()

H_0.25 <- readRDS("output/RDS/ZEB2expressingGenes_eachCluster_0.25_Human.rds")
C_0.25 <- readRDS("output/RDS/ZEB2expressingGenes_eachCluster_0.25_Chimp.rds")
H_0.5 <- readRDS("output/RDS/ZEB2expressingGenes_eachCluster_0.5_Human.rds")
C_0.5 <- readRDS("output/RDS/ZEB2expressingGenes_eachCluster_0.5_Chimp.rds")

ZEB2exp <- c(H_0.5, C_0.5) %>% as_tibble() %>% 
  add_column(species = rep(c("Human", "Chimp"), each =11 )) %>%
  add_column(CT = rep(names(table(Idents(integrated))), 2))
ZEB2exp$CT <- factor(ZEB2exp$CT, levels = names(table(Idents(integrated))))
ZEB2exp$species <- factor(ZEB2exp$species, levels = c("Human", "Chimp"))

#Number of ZEB2 expressing samples
ggplot(data=ZEB2exp, aes(x=CT, y=value, color=species)) +
  geom_point(aes(color=species), size=3, stat="identity") +
  scale_y_continuous(limits=c(0, 9400), breaks=c(0, 500, 1000, 2000, 5000, 9000))+
  scale_y_break(c(2000, 3000, 4000, 9000)) +
  theme(axis.title = element_blank()) +
  labs(title="Number of ZEB2 expressing samples in each cluster", x="", y="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_color_manual(values=c("#8dd35fff", "#ff80b2ff"))

#barplot showing the percentage of correlated genes against the total number of ZEB2 expressing genes
#if the group_by doesn't work, unload "plyr"

#Nrcorr_per <- nr_corr_sort
#order the levels
#Nrcorr_per$species <- factor(Nrcorr_per$species, levels =c("human", "chimpanzee"))
#Nrcorr_per$CT <- factor(Nrcorr_per$CT, levels = names(table(Idents(integrated))))
#Nrcorr_per <- Nrcorr_per %>% arrange(species, CT) %>% mutate(n_ZEBsamples = rep(ZEB2exp$value, each = 2))

#H_0.5_ALL <- readRDS("output/RDS/ExpressedGenes_eachCluster_0.5_Human.rds")
#C_0.5_ALL <- readRDS("output/RDS/ExpressedGenes_eachCluster_0.5_Chimp.rds")
#This was a wrong approach!!!

#number of cells in each cluster
table(integrated@meta.data$seurat_clusters)

#number of cells in each cluster*species
table(integrated@meta.data[which(stringr::str_detect(rownames(integrated@meta.data), "Human")),24])
table(integrated@meta.data[which(stringr::str_detect(rownames(integrated@meta.data), "Chimp")),24])

ExpPattern <- ZEB2exp %>% mutate(totalcells = c(table(integrated@meta.data[which(stringr::str_detect(rownames(integrated@meta.data), "Human")),24]),
                                                table(integrated@meta.data[which(stringr::str_detect(rownames(integrated@meta.data), "Chimp")),24])))
ExpPattern <- ExpPattern %>% mutate(per = value/totalcells*100)

#proportion of ZEB2 expressing cells 
ggplot(data=filter(ExpPattern), aes(x=CT, y=per, fill=species)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= per, label=round(per)), vjust=1.3, color="white", size=2.5)+
  scale_y_continuous(limits=c(0, 100))+
  facet_wrap( ~ species) +
  theme(axis.title = element_blank()) +
  labs(title="Proportion of ZEB2 expressing cells", x="", y="percentage (%)") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_manual(values=c("#8dd35fff", "#ff80b2ff")) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1)) 


#############################################################################
#a venn diagram showing the overlapping genes between human & chimp without considering CT
#############################################################################
#make a function to collect the gene names of every significant cor.test
GetGeneName.pos <-function(x){ 
  n <- c()
  for(i in 1:ncol(x))
  {
    n <- c(n, rownames(x)[which(sign(x[,i])==1)])
  }
  return(n)
}

GetGeneName.neg <-function(x){ 
  n <- c()
  for(i in 1:ncol(x))
  {
    n <- c(n, rownames(x)[which(sign(x[,i])==-1)])
  }
  return(n)
}


length(GetGeneName.pos(cor_H)) == sum(nr_corr$n[1:11]) #467 - not unique!!!
length(GetGeneName.pos(cor_C)) == sum(nr_corr$n[23:33]) #1150 - not unique!!!

length(GetGeneName.neg(cor_H)) == sum(nr_corr$n[12:22]) #602 - not unique!!!
length(GetGeneName.neg(cor_C)) == sum(nr_corr$n[34:44]) #402 - not unique!!!

#merge pos and neg into a list
GeneNames.human <- c(GetGeneName.pos(cor_H), GetGeneName.neg(cor_H))
GeneNames.chimp <- c(GetGeneName.pos(cor_C), GetGeneName.neg(cor_C))

#ven diagram library
library(ggvenn)
length(intersect(GeneNames.human, GeneNames.chimp)) #313


#Create the input data for venn diagram
overlaps.hvsc.pn <- list(
  pos.human = GetGeneName.pos(cor_H),
  neg.human = GetGeneName.neg(cor_H),
  neg.chimp = GetGeneName.neg(cor_C),
  pos.chimp = GetGeneName.pos(cor_C)
)

overlaps.hvsc <- list(
  human = c(GetGeneName.pos(cor_H),GetGeneName.neg(cor_H)),
  chimp = c(GetGeneName.pos(cor_C),GetGeneName.neg(cor_C))
)

CTList_H <- lapply(1:ncol(cor_H), function(x){x <- rownames(cor_H)[which(!is.na(cor_H[,x]))]})
names(CTList_H) <- paste0(names(table(Idents(integrated))))
CTList_C <- lapply(1:ncol(cor_C), function(x){x <- rownames(cor_C)[which(!is.na(cor_C[,x]))]})
names(CTList_C) <- paste0(names(table(Idents(integrated))))

#overlapping genes with top 3 highest number of correlated genes in human
Reduce(intersect, list(CTList_H[[1]],CTList_H[[2]], CTList_H[[8]]))
#overlapping genes with top 2 highest number of correlated genes in chimp
Reduce(intersect, list(CTList_C[[1]],CTList_H[[10]]))

intersect(CTList_H[[1]],CTList_H[[2]])

#I manually created plotList1 to 11 :P
plotList1 <- lapply(
  1:ncol(cor_H),
  function(CT) {
    x <- ggvenn(
      CTList_H[c(1,CT)],
      fill_color = c(unname(my_cols)[1], unname(my_cols)[CT]),
      show_percentage = FALSE
    )
    x
  }
)

plotList1 <- lapply(
  1:ncol(cor_C),
  function(CT) {
    x <- ggvenn(
      CTList_C[c(1,CT)],
      fill_color = c(unname(my_cols)[1], unname(my_cols)[CT]),
      show_percentage = FALSE
    )
    x
  }
)

allplots <- ggarrange(plotlist=c(plotList1, plotList2, plotList3, plotList4, plotList5,
                                 plotList6, plotList7, plotList8, plotList9,plotList10,plotList11), 
                      ncol = 11, nrow = 11)

allplots

#color setting
#neg col: #66c2a4 - #648077
#pos.col: #fc8d62 - #927063

#compute the venn diagram
ggvenn(
  overlaps.hvsc.pn, 
  fill_color = c("#fc8d62", "#66c2a4", "#648077", "#927063"),
  stroke_size = 0.5, set_name_size = 4
) + labs(title="Number of correlated genes (+/-)") +
  theme(plot.title = element_text(hjust = 0.5))

ggvenn(
  overlaps.hvsc, 
  fill_color = c("#8dd35fff", "#ff80b2ff"),
  stroke_size = 0.5, set_name_size = 4
) + labs(title="Number of correlated genes") +
  theme(plot.title = element_text(hjust = 0.5))

length(setdiff(c(GetGeneName.pos(cor_H),GetGeneName.neg(cor_H)), c(GetGeneName.pos(cor_C),GetGeneName.neg(cor_C)))) #502 h-specific - already unique
length(setdiff(c(GetGeneName.pos(cor_C),GetGeneName.neg(cor_C)), c(GetGeneName.pos(cor_H),GetGeneName.neg(cor_H)))) #1141 c-specific - already unique
length(intersect(c(GetGeneName.pos(cor_H),GetGeneName.neg(cor_H)), c(GetGeneName.pos(cor_C),GetGeneName.neg(cor_C)))) #313 overlaps - already unique

ALLoverlaps <- intersect(c(GetGeneName.pos(cor_H),GetGeneName.neg(cor_H)), c(GetGeneName.pos(cor_C),GetGeneName.neg(cor_C)))
write.table(ALLoverlaps, file = paste0(path_out, "OverlappingCorrelatedGenes_allgenes.tsv"), 
            sep = "\t", row.names = F, col.names = T, quote = F)


#############################################################################
#comparison of correlated genes of human and chimp in cluster 0, 1, 3, 7, 9
#############################################################################
CTList <- c("0", "1", "3", "7", "9")
plotList_CT <- lapply(
  1:length(CTList),
  function(CT) {
    input <- list(CTList_H[[CTList[CT]]], CTList_C[[CTList[CT]]])
    names(input) <- c("human", "chimp")
    x <- ggvenn(
      input,
      fill_color = c("#8dd35fff", "#ff80b2ff"),
      show_percentage = FALSE
    ) + labs(title=paste0("cluster", CTList[CT])) +
      theme(plot.title = element_text(hjust = 0.5))
    x
  }
)

allplots_CT <- ggarrange(plotlist=plotList_CT, 
                      ncol = 5, nrow = 1)

allplots_CT

length(intersect(CTList_H[[1]], CTList_C[[1]])) #313 overlaps - already unique

#export the gene list of each overlapping genes per clusters
for(i in 1:length(CTList)){
  overlap <- intersect(CTList_H[[CTList[i]]], CTList_C[[CTList[i]]])
  write.table(overlap, file = paste0(path_out, 
                                     "OverlappingCorrelatedGenes_CT_", CTList[i],
                                     "_n", length(intersect(CTList_H[[CTList[i]]], CTList_C[[CTList[i]]])), ".tsv"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
}

#export the gene list of each h-specific genes per clusters
for(i in 1:length(CTList)){
  overlap <- setdiff(CTList_H[[CTList[i]]], CTList_C[[CTList[i]]])
  write.table(overlap, file = paste0(path_out, 
                                     "Hspe_correlatedgenes_CT_", CTList[i], 
                                    "_n",  length(setdiff(CTList_H[[CTList[i]]], CTList_C[[CTList[i]]])),".tsv"), 
              sep = "\t", row.names = F, col.names = T, quote = F)
}

setdiff(CTList_H[[CTList[4]]], CTList_C[[CTList[4]]])

#############################################################################
#human-specific genes
#############################################################################
hspec <- setdiff(c(GetGeneName.pos(cor_H),GetGeneName.neg(cor_H)), c(GetGeneName.pos(cor_C),GetGeneName.neg(cor_C)))
#table of human-specific genes with correlation coefficient
tmp <- cor_H[which(rownames(cor_H) %in% hspec),]

nr_corr_hspec <- data.frame(matrix(, nrow=c(ncol(cor_H))*2, ncol=3))
colnames(nr_corr_hspec) <- c("CT", "correlation","n")
nr_corr_hspec[,1] <- c(rep(colnames(cor_H), 2))
nr_corr_hspec[,2] <- c(rep("pos", ncol(cor_H)), rep("neg",ncol(cor_H)))

for (i in 1:nrow(nr_corr_hspec)){
    if (i < ncol(cor_H)+1){
      nr_corr_hspec[i,3] <- as.numeric(length(rownames(tmp)[which(sign(tmp[,i]) == 1)]))
    } else {
      nr_corr_hspec[i,3] <- as.numeric(length(rownames(tmp)[which(sign(tmp[,i-11])==-1)]))
    }
}

#change the nr of corr into numeric
nr_corr_hspec$n <- as.numeric(nr_corr_hspec$n)

library(plyr)
# Sort by timpoint
nr_corr_hspec_sort <- arrange(nr_corr_hspec, CT)
# Calculate the cumulative sum of n for each correlation
nr_corr_hspec_sum <- ddply(nr_corr_hspec_sort, "CT",
                     transform, label_ypos=cumsum(n))

#calculate the log values of counts
nr_corr_hspec_sum$log_n <- log(nr_corr_hspec_sum$n)
#change the infinite values to zero
nr_corr_hspec_sum$log_n[which(!is.finite(nr_corr_hspec_sum$log_n))] <- 0
#calculate the log values of cumulative sums
nr_corr_hspec_sum$log_label_ypos <- log(nr_corr_hspec_sum$label_ypos)

#order the levels
nr_corr_hspec_sum$CT <- factor(nr_corr_hspec_sum$CT, levels =colnames(cor_H))

#barplot with number of correlated genes
ggplot(data=filter(nr_corr_hspec_sum), aes(x=CT, y=n, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.3, color="white", size=2.5)+
  scale_y_continuous(limits=c(0, 350), breaks = c(0, 100, 200, 300)) +
  theme(axis.title = element_blank()) +
  labs(title="Number of correlated genes with ZEB2 expression (human-specific)", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2") + 
  scale_y_break(c(50, 75, 150, 300), scales = 0.1) +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))

#distribution of human specific correlated genes in each cluster
ggplot(data=filter(nr_corr_hspec_sum), aes(x = correlation, y=n, fill=CT)) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.title = element_blank()) +
  labs(title="Proportion of correlated genes in each cluster (human-specific)", x="", y="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))


#human-specific genes per cluster
GetGeneName.pos.ct <- function(x, y){ 
  n <- c()
  n <- c(n, rownames(x)[which(sign(x[,y])==1)])
  return(n)
}

GetGeneName.neg.ct <- function(x, y){ 
  n <- c()
  n <- c(n, rownames(x)[which(sign(x[,y])==-1)])
  return(n)
}
  
for (i in 1:ncol(tmp)) {
  x <- data.frame("gene_symbol" = GetGeneName.pos.ct(tmp, colnames(tmp)[i]))
  write.table(x, file = paste(path_out,'Hspe_correlatedgenes_pos_CT_',colnames(tmp)[i],'_n',length(GetGeneName.pos(tmp, colnames(tmp)[i])),'.tsv', sep = ""), row.names = F, col.names = T, quote = F)
  x <- data.frame("gene_symbol" = GetGeneName.neg.ct(tmp, colnames(tmp)[i]))
  write.table(x, file = paste(path_out,'Hspe_correlatedgenes_neg_CT_',colnames(tmp)[i],'_n',length(GetGeneName.neg(tmp, colnames(tmp)[i])),'.tsv', sep = ""), row.names = F, col.names = T, quote = F)
}


#############################################################################
#GO based on h & c specific without considering CT
#############################################################################
#library calling
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#GOs for three groups - overlaps, human-specific, chimp-specific

#overlapping genes
overlaps <- as.data.frame(intersect(c(GetGeneName.pos(cor_H),GetGeneName.neg(cor_H)), c(GetGeneName.pos(cor_C),GetGeneName.neg(cor_C))))
colnames(overlaps)[1] <- "ID"
overlaps$group <- "overlaps"

#human-specific
cluster_h_sp <- as.data.frame(setdiff(c(GetGeneName.pos(cor_H),GetGeneName.neg(cor_H)), c(GetGeneName.pos(cor_C),GetGeneName.neg(cor_C))))
colnames(cluster_h_sp)[1] <- "ID"
cluster_h_sp$group <- "human-specific"

#gorilla-specific
cluster_c_sp <- as.data.frame(setdiff(c(GetGeneName.pos(cor_C),GetGeneName.neg(cor_C)), c(GetGeneName.pos(cor_H),GetGeneName.neg(cor_H))))
colnames(cluster_c_sp)[1] <- "ID"
cluster_c_sp$group <- "chimp-specific"


cluster_over.go <- enrichGO(gene=overlaps[,1],
                            keyType = "SYMBOL", 
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 1)

cluster_h_sp.go <- enrichGO(gene=cluster_h_sp[,1],
                            keyType = "SYMBOL", 
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "none",
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 1)

cluster_c_sp.go <- enrichGO(gene=cluster_c_sp[,1],
                            keyType = "SYMBOL", 
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "none",
                            pvalueCutoff = 0.05, 
                            qvalueCutoff = 1)

dotplot(cluster_over.go, showCategory=15) + ggtitle("overlapping cluster")
dotplot(cluster_h_sp.go, showCategory=15) + ggtitle("human-specific cluster")
dotplot(cluster_c_sp.go, showCategory=15) + ggtitle("chimp-specific cluster")

#dotplot
library(scales)
p1 <- dotplot(cluster_over.go, showCategory=15) + ggtitle("overlapping cluster")
p1 <- p1 + scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none') 
p2 <- dotplot(cluster_h_sp.go, showCategory=15) + ggtitle("human-specific cluster")
p2 <- p2+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
p3 <- dotplot(cluster_c_sp.go, showCategory=15) + ggtitle("chimp-specific cluster")
p3 <- p3+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')

cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3])

#emapplot
cluster_merg <- rbind(overlaps, cluster_h_sp, cluster_c_sp)

merg_res <- compareCluster(ID~group, data=cluster_merg, fun="enrichGO", OrgDb = org.Hs.eg.db,keyType = "SYMBOL", pvalueCutoff=0.05)
merg_res <- pairwise_termsim(merg_res)
emapplot(merg_res) 

#identity of genes specific to these three clusters within the GO "cadherin binding" according to merg_res emapplot result.
#are they really separate genes, or they are for instance paralogs that were separated by our liftover to human genome.
cad_genes <- merg_res@compareClusterResult$geneID[which(merg_res@compareClusterResult$Description == "cadherin binding")]
names(cad_genes) <- merg_res@compareClusterResult$group[which(merg_res@compareClusterResult$Description == "cadherin binding")]
cad_genes <- strsplit(cad_genes, "/")

for (i in 1:3){
  write.csv(cad_genes[i], file=paste('cadherin_genes_emapplot_',names(cad_genes)[i],'.csv', sep = ""), row.names=F)
} 

#another aspect from later future on January 10th 2023
#species-specific pos vs neg GO

h_pos <- data.frame(geneID =setdiff(setdiff(setdiff(GetGeneName.pos(cor_H), GetGeneName.neg(cor_H)),  GetGeneName.neg(cor_C)), GetGeneName.pos(cor_C)),
                    coe = "(+)", species = "Hu")
h_neg <- data.frame(geneID =setdiff(setdiff(setdiff(GetGeneName.neg(cor_H), GetGeneName.pos(cor_H)),  GetGeneName.neg(cor_C)), GetGeneName.pos(cor_C)),
                    coe = "(-)", species = "Hu")
c_pos <- data.frame(geneID =setdiff(setdiff(setdiff(GetGeneName.pos(cor_C), GetGeneName.neg(cor_C)),  GetGeneName.neg(cor_H)), GetGeneName.pos(cor_H)),
                    coe = "(+)", species = "Ch")
c_neg <- data.frame(geneID =setdiff(setdiff(setdiff(GetGeneName.neg(cor_C), GetGeneName.pos(cor_C)),  GetGeneName.neg(cor_H)), GetGeneName.pos(cor_H)),
                    coe = "(-)", species = "Ch")
hc_pn_GO <- rbind(h_pos, h_neg, c_pos, c_neg)

#all genes either positively or negatively correlated with ZEB2 
write.table(GetGeneName.pos(cor_H), file = paste(path_out,'correlatedgenes_H_pos_n',length(GetGeneName.pos(cor_H)), '.txt', sep = ""), row.names = F, col.names = T, quote = F)
write.table(GetGeneName.neg(cor_H), file = paste(path_out,'correlatedgenes_H_neg_n',length(GetGeneName.neg(cor_H)),'.txt', sep = ""), row.names = F, col.names = T, quote = F)
write.table(GetGeneName.pos(cor_C), file = paste(path_out,'correlatedgenes_C_pos_n',length(GetGeneName.pos(cor_C)),'.txt', sep = ""), row.names = F, col.names = T, quote = F)
write.table(GetGeneName.neg(cor_C), file = paste(path_out,'correlatedgenes_C_neg_n',length(GetGeneName.neg(cor_C)),'.txt', sep = ""), row.names = F, col.names = T, quote = F)

#human-specifically positively or negatively correlated genes
write.table(h_pos$geneID, file = paste(path_out,'correlatedgenes_H_pos_specific_n',nrow(h_pos), '.txt', sep = ""), row.names = F, col.names = F, quote = F)
write.table(h_neg$geneID, file = paste(path_out,'correlatedgenes_H_neg_specific_n',nrow(h_neg), '.txt', sep = ""), row.names = F, col.names = F, quote = F)


compare.formular_HC <- compareCluster(geneID~species+coe, 
                                             data =hc_pn_GO,
                                             pvalueCutoff = 0.05,
                                             fun = "enrichGO",
                                             OrgDb="org.Hs.eg.db", 
                                             keyType = "ENSEMBL", 
                                             ont = "BP")
#No enrichment found in any of gene cluster, please check your input...
compare.formular_HC <- compareCluster(geneID~species, 
                                      data =hc_pn_GO,
                                      pvalueCutoff = 0.05,
                                      fun = "enrichGO",
                                      OrgDb="org.Hs.eg.db", 
                                      keyType = "ENSEMBL", 
                                      ont = "BP")
#No enrichment found in any of gene cluster, please check your input...

#############################################################################
#Plotting the number of correlated genes with ZEB2 in three overlapping CT
#############################################################################
#overlapping CT: cortical neurons 1, cortical neurons 2, mesenchymal-like cells

nr_corr_CT3 <- rbind(nr_corr[which(str_detect(nr_corr[,1], "cortical neurons 1")),],
                     nr_corr[which(str_detect(nr_corr[,1], "cortical neurons 2")),],
                     nr_corr[which(str_detect(nr_corr[,1], "mesenchymal-like cells")),])

#change the nr of corr into numeric
nr_corr_CT3$n <- as.numeric(nr_corr_CT3$n)

# Sort by timpoint and species
nr_corr_CT3_sort <- arrange(nr_corr_CT3, CT, species)
# Calculate the cumulative sum of n for each correlation
nr_corr_CT3_sum <- ddply(nr_corr_CT3_sort, c("CT", "species"),
                     transform, label_ypos=cumsum(n))

ggplot(data=nr_corr_CT3_sum, aes(x=CT, y=n, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.3, color="white", size=2.5)+
  scale_y_continuous(limits=c(0, 300), breaks=c(0, 100, 200, 300))+
  facet_wrap( ~ species) +
  theme(axis.title = element_blank()) +
  labs(title="Number of correlated genes with ZEB2 expression \nin cell types expressing ZEB2 in both species", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 18))+
  scale_fill_brewer(palette = "Set2") + 
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 60, hjust = 1))

#Size of ggplot above: 950x530

#############################################################################
#a venn diagram showing the overlapping genes between human & chimp with considering CT, aka in three CTs
#############################################################################
#make a function to collect the gene names of every significant cor.test
GetGeneName.pos <-function(x, y){ 
  n <- c()
  n <- c(n, rownames(x)[which(sign(x[,y])==1)])
  return(n)
}

GetGeneName.neg <-function(x, y){ 
  n <- c()
  n <- c(n, rownames(x)[which(sign(x[,y])==-1)])
  return(n)
}

#Create the input data for venn diagram
overlaps.hvsc.pn_CN1 <- list(
  human.pos =GetGeneName.pos(cor_H, "cortical neurons 1"),
  human.neg =GetGeneName.neg(cor_H, "cortical neurons 1"),
  chimp.neg =GetGeneName.neg(cor_C, "cortical neurons 1"),
  chimp.pos =GetGeneName.pos(cor_C, "cortical neurons 1")
)

overlaps.hvsc.pn_CN2 <- list(
  human.pos =GetGeneName.pos(cor_H, "cortical neurons 2"),
  human.neg =GetGeneName.neg(cor_H, "cortical neurons 2"),
  chimp.neg =GetGeneName.neg(cor_C, "cortical neurons 2"),
  chimp.pos =GetGeneName.pos(cor_C, "cortical neurons 2")
)

overlaps.hvsc.pn_MLC <- list(
  human.pos =GetGeneName.pos(cor_H, "mesenchymal-like cells"),
  human.neg =GetGeneName.neg(cor_H, "mesenchymal-like cells"),
  chimp.neg =GetGeneName.neg(cor_C, "mesenchymal-like cells"),
  chimp.pos =GetGeneName.pos(cor_C, "mesenchymal-like cells")
)


#color setting
#neg col: #66c2a4 - #648077
#pos.col: #fc8d62 - #927063

#compute the venn diagram
v1 <- ggvenn(
  overlaps.hvsc.pn_CN1, 
  fill_color = c("#fc8d62", "#66c2a4", "#648077", "#927063"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE) +
  ggtitle("cortical neurons 1") + 
  theme(plot.title = element_text(size=18,  face="bold",hjust = 0.5))

v2 <- ggvenn(
  overlaps.hvsc.pn_CN2, 
  fill_color = c("#fc8d62", "#66c2a4", "#648077", "#927063"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE)+
  ggtitle("cortical neurons 2") + 
  theme(plot.title = element_text(size=18,  face="bold",hjust = 0.5))

v3 <- ggvenn(
  overlaps.hvsc.pn_MLC, 
  fill_color = c("#fc8d62", "#66c2a4", "#648077", "#927063"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE)+
  ggtitle("mesenchymal-like cells") + 
  theme(plot.title = element_text(size=18,  face="bold",hjust = 0.5))


cowplot::plot_grid(v1, v2, v3, ncol=3, labels=LETTERS[1:3])

#only three overlaps 
#1) cortical neurons 1 (human.neg & chimp.neg)
#2) cortical neurons 1 (human.pos & chimp.pos)
#3) cortical neurons 2 (human.pos & chimp.pos)
#4) mesenchymal-like cells (human.neg & chimp.pos)
#5) mesenchymal-like cells (human.pos & chimp.neg)
#6) mesenchymal-like cells (human.pos & chimp.pos)

overlap_1 <- intersect(GetGeneName.neg(cor_H, "cortical neurons 1"), GetGeneName.neg(cor_C, "cortical neurons 1"))
overlap_2 <- intersect(GetGeneName.pos(cor_H, "cortical neurons 1"), GetGeneName.pos(cor_C, "cortical neurons 1"))

overlap_3 <- intersect(GetGeneName.neg(cor_H, "mesenchymal-like cells"), GetGeneName.pos(cor_C, "mesenchymal-like cells"))
overlap_4 <- intersect(GetGeneName.pos(cor_H, "mesenchymal-like cells"), GetGeneName.neg(cor_C, "mesenchymal-like cells"))

#export the gene lists
ct_names <- c("cortical neurons 1", "cortical neurons 2", "mesenchymal-like cells")
#all genes either positively or negatively correlated with ZEB2 in three cell types 
#human
for (i in 1:length(ct_names)) {
  tmp <- data.frame("gene symbol" = GetGeneName.pos(cor_H, ct_names[i]))
  write.table(tmp, file = paste(path_out,'correlatedgenes_pos_H_',ct_names[i],'_n',length(GetGeneName.pos(cor_H, ct_names[i])),'.txt', sep = ""), row.names = F, col.names = F, quote = F)
  tmp <- data.frame("gene symbol" = GetGeneName.neg(cor_H, ct_names[i]))
  write.table(tmp, file = paste(path_out,'correlatedgenes_neg_H_',ct_names[i],'_n',length(GetGeneName.neg(cor_H, ct_names[i])),'.txt', sep = ""), row.names = F, col.names = F, quote = F)
}

#chimp
for (i in 1:length(ct_names)) {
  tmp <- data.frame("gene symbol" = GetGeneName.pos(cor_C, ct_names[i]))
  write.table(tmp, file = paste(path_out,'correlatedgenes_pos_C_',ct_names[i],'_n',length(GetGeneName.pos(cor_C, ct_names[i])),'.txt', sep = ""), row.names = F, col.names = F, quote = F)
  tmp <- data.frame("gene symbol" = GetGeneName.neg(cor_C, ct_names[i]))
  write.table(tmp, file = paste(path_out,'correlatedgenes_neg_C_',ct_names[i],'_n',length(GetGeneName.neg(cor_C, ct_names[i])),'.txt', sep = ""), row.names = F, col.names = F, quote = F)
}

#human-specifically positively or negatively correlated genes
write.table(setdiff(GetGeneName.pos(cor_H, ct_names[1]), GetGeneName.pos(cor_C, ct_names[1])), 
            file = paste(path_out,'correlatedgenes_H_pos_specific_CT1_n',length(setdiff(GetGeneName.pos(cor_H, ct_names[1]), GetGeneName.pos(cor_C, ct_names[1]))), '.txt', sep = ""), 
            row.names = F, col.names = F, quote = F)
write.table(setdiff(GetGeneName.neg(cor_H, ct_names[1]), GetGeneName.neg(cor_C, ct_names[1])), 
            file = paste(path_out,'correlatedgenes_H_neg_specific_CT1_n',length(setdiff(GetGeneName.neg(cor_H, ct_names[1]), GetGeneName.neg(cor_C, ct_names[1]))), '.txt', sep = ""), 
            row.names = F, col.names = F, quote = F)

write.table(setdiff(GetGeneName.pos(cor_H, ct_names[3]), GetGeneName.neg(cor_C, ct_names[3])), 
            file = paste(path_out,'correlatedgenes_H_pos_specific_MLC_n',length(setdiff(GetGeneName.pos(cor_H, ct_names[3]), GetGeneName.neg(cor_C, ct_names[3]))), '.txt', sep = ""), 
            row.names = F, col.names = F, quote = F)
write.table(setdiff(GetGeneName.neg(cor_H, ct_names[3]), GetGeneName.pos(cor_C, ct_names[3])), 
            file = paste(path_out,'correlatedgenes_H_neg_specific_MLC_n',length(setdiff(GetGeneName.neg(cor_H, ct_names[3]), GetGeneName.pos(cor_C, ct_names[3]))), '.txt', sep = ""), 
            row.names = F, col.names = F, quote = F)


setdiff(GetGeneName.pos(cor_H, ct_names[1]), GetGeneName.pos(cor_C, ct_names[1]))

#############################################################################
#GO based on h & c specific without in three CTs
#############################################################################
#GOs for three groups - human.pos-specific, human.neg-specific, chimp.pos-specific, chimp.neg-specific

#cortical neurons 1
#human.pos-specific
human.pos.CT1 <- as.data.frame(setdiff(c(GetGeneName.pos(cor_H, "cortical neurons 1")), c(GetGeneName.pos(cor_C, "cortical neurons 1"))))
colnames(human.pos.CT1)[1] <- "ID"
human.pos.CT1$group <- "human.pos"

#human.neg-specific
human.neg.CT1 <- as.data.frame(setdiff(c(GetGeneName.neg(cor_H, "cortical neurons 1")), c(GetGeneName.neg(cor_C, "cortical neurons 1"))))
colnames(human.neg.CT1)[1] <- "ID"
human.neg.CT1$group <- "human.neg"

#chimp.pos-specific
chimp.pos.CT1 <- as.data.frame(setdiff(c(GetGeneName.pos(cor_C, "cortical neurons 1")), c(GetGeneName.pos(cor_H, "cortical neurons 1"))))
colnames(chimp.pos.CT1)[1] <- "ID"
chimp.pos.CT1$group <- "chimp.pos"

#chimp.pos-specific
chimp.neg.CT1 <- as.data.frame(setdiff(c(GetGeneName.neg(cor_C, "cortical neurons 1")), c(GetGeneName.neg(cor_H, "cortical neurons 1"))))
colnames(chimp.neg.CT1)[1] <- "ID"
chimp.neg.CT1$group <- "chimp.neg"


human.pos.CT1.go <- enrichGO(gene=human.pos.CT1[,1],
                            keyType = "SYMBOL", 
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01, 
                            qvalueCutoff = 0.05)
human.pos.CT1.go@result$geneID

human.neg.CT1.go <- enrichGO(gene=human.neg.CT1[,1],
                            keyType = "SYMBOL", 
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01, 
                            qvalueCutoff = 0.05)

human.neg.CT1.go@result$geneID[1:15]

chimp.pos.CT1.go <- enrichGO(gene=chimp.pos.CT1[,1],
                            keyType = "SYMBOL", 
                            OrgDb = org.Hs.eg.db, 
                            ont = "BP", 
                            pAdjustMethod = "BH",
                            pvalueCutoff = 0.01, 
                            qvalueCutoff = 0.05)

chimp.neg.CT1.go <- enrichGO(gene=chimp.neg.CT1[,1],
                             keyType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.01, 
                             qvalueCutoff = 0.05)


d <- godata('org.Hs.eg.db', ont="BP")

human.pos.CT1.go_test <- pairwise_termsim(human.pos.CT1.go, method = "Wang", semData = d)
tt <- treeplot(human.pos.CT1.go_test)

human.neg.CT1.go_test <- pairwise_termsim(human.neg.CT1.go, method = "Wang", semData = d)
tt2 <- treeplot(human.neg.CT1.go_test)

barplot(human.pos.CT1.go, showCategory=15) + ggtitle("human.pos-specific cluster")
barplot(human.neg.CT1.go, showCategory=15) + ggtitle("human.neg-specific cluster")
barplot(chimp.pos.CT1.go, showCategory=15) + ggtitle("chimp.pos-specific cluster") # 0 enrichment
barplot(chimp.neg.CT1.go, showCategory=15) + ggtitle("chimp.neg-specific cluster") # 0 enrichment

#dotplot
library(scales)
ct11 <- dotplot(human.pos.CT1.go, showCategory=15) + ggtitle("human.pos-specific cluster")
ct11 <- ct11 + scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none') 
ct12 <- dotplot(human.neg.CT1.go, showCategory=15) + ggtitle("human.neg-specific cluster")
ct12 <- ct12+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
ct13 <- dotplot(chimp.pos.CT1.go, showCategory=15) + ggtitle("chimp.pos-specific cluster")
ct13 <- ct13+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
ct14 <- dotplot(chimp.neg.CT1.go, showCategory=15) + ggtitle("chimp.neg-specific cluster")
ct14 <- ct14+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')


cowplot::plot_grid(ct11, ct12, ct13, ct14, ncol=2, labels=LETTERS[1:4])

#emapplot
ct1_merg <- rbind(human.pos.CT1, human.neg.CT1, chimp.pos.CT1, chimp.neg.CT1)
write.table(ct1_merg, file = "correlatedgenes_species-specific_CT1.txt", row.names = F, quote = F)

ct1_merg_res <- compareCluster(ID~group, data=ct1_merg, fun="enrichGO", OrgDb = org.Hs.eg.db,keyType = "SYMBOL", pvalueCutoff=0.05)
ct1_merg_res <- pairwise_termsim(ct1_merg_res)
emapplot(ct1_merg_res) 

#cortical neurons 2
#human.pos-specific
human.pos.CT2 <- as.data.frame(setdiff(c(GetGeneName.pos(cor_H, "cortical neurons 2")), c(GetGeneName.pos(cor_C, "cortical neurons 2"))))
colnames(human.pos.CT2)[1] <- "ID"
human.pos.CT2$group <- "human.pos"

#human.neg-specific
human.neg.CT2 <- as.data.frame(setdiff(c(GetGeneName.neg(cor_H, "cortical neurons 2")), c(GetGeneName.neg(cor_C, "cortical neurons 2"))))
colnames(human.neg.CT2)[1] <- "ID"
human.neg.CT2$group <- "human.neg"

#chimp.pos-specific
chimp.pos.CT2 <- as.data.frame(setdiff(c(GetGeneName.pos(cor_C, "cortical neurons 2")), c(GetGeneName.pos(cor_H, "cortical neurons 2"))))
colnames(chimp.pos.CT2)[1] <- "ID"
chimp.pos.CT2$group <- "chimp.pos"

#chimp.pos-specific
chimp.neg.CT2 <- as.data.frame(setdiff(c(GetGeneName.neg(cor_C, "cortical neurons 2")), c(GetGeneName.neg(cor_H, "cortical neurons 2"))))
colnames(chimp.neg.CT2)[1] <- "ID"
chimp.neg.CT2$group <- "chimp.neg"


human.pos.CT2.go <- enrichGO(gene=human.pos.CT2[,1],
                             keyType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 1)

human.neg.CT2.go <- enrichGO(gene=human.neg.CT2[,1],
                             keyType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "none",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 1)

chimp.pos.CT2.go <- enrichGO(gene=chimp.pos.CT2[,1],
                             keyType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "none",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 1)

chimp.neg.CT2.go <- enrichGO(gene=chimp.neg.CT2[,1],
                             keyType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "none",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 1)

dotplot(human.pos.CT2.go, showCategory=15) + ggtitle("human.pos-specific cluster")
dotplot(human.neg.CT2.go, showCategory=15) + ggtitle("human.neg-specific cluster")
dotplot(chimp.pos.CT2.go, showCategory=15) + ggtitle("chimp.pos-specific cluster")
dotplot(chimp.neg.CT2.go, showCategory=15) + ggtitle("chimp.neg-specific cluster")

#dotplot
CT21 <- dotplot(human.pos.CT2.go, showCategory=15) + ggtitle("human.pos-specific cluster")
CT21 <- CT21 + scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none') 
CT22 <- dotplot(human.neg.CT2.go, showCategory=15) + ggtitle("human.neg-specific cluster")
CT22 <- CT22+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
CT23 <- dotplot(chimp.pos.CT2.go, showCategory=15) + ggtitle("chimp.pos-specific cluster")
CT23 <- CT23+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
CT24 <- dotplot(chimp.neg.CT2.go, showCategory=15) + ggtitle("chimp.neg-specific cluster")
CT24 <- CT24+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')


cowplot::plot_grid(CT21, CT22, CT23, CT24, ncol=2, labels=LETTERS[1:4])

#emapplot
CT2_merg <- rbind(human.pos.CT2, human.neg.CT2, chimp.pos.CT2, chimp.neg.CT2)
write.table(CT2_merg, file = "correlatedgenes_species-specific_CT2.txt", row.names = F, quote = F)

CT2_merg_res <- compareCluster(ID~group, data=CT2_merg, fun="enrichGO", OrgDb = org.Hs.eg.db,keyType = "SYMBOL", pvalueCutoff=0.05)
CT2_merg_res <- pairwise_termsim(CT2_merg_res)
emapplot(CT2_merg_res) 


#mesenchymal-like cells
#human.pos-specific
human.pos.MCL <- as.data.frame(setdiff(c(GetGeneName.pos(cor_H, "mesenchymal-like cells")), c(GetGeneName.pos(cor_C, "mesenchymal-like cells"))))
colnames(human.pos.MCL)[1] <- "ID"
human.pos.MCL$group <- "human.pos"

#human.neg-specific
human.neg.MCL <- as.data.frame(setdiff(c(GetGeneName.neg(cor_H, "mesenchymal-like cells")), c(GetGeneName.neg(cor_C, "mesenchymal-like cells"))))
colnames(human.neg.MCL)[1] <- "ID"
human.neg.MCL$group <- "human.neg"

#chimp.pos-specific
chimp.pos.MCL <- as.data.frame(setdiff(c(GetGeneName.pos(cor_C, "mesenchymal-like cells")), c(GetGeneName.pos(cor_H, "mesenchymal-like cells"))))
colnames(chimp.pos.MCL)[1] <- "ID"
chimp.pos.MCL$group <- "chimp.pos"

#chimp.pos-specific
chimp.neg.MCL <- as.data.frame(setdiff(c(GetGeneName.neg(cor_C, "mesenchymal-like cells")), c(GetGeneName.neg(cor_H, "mesenchymal-like cells"))))
colnames(chimp.neg.MCL)[1] <- "ID"
chimp.neg.MCL$group <- "chimp.neg"


human.pos.MCL.go <- enrichGO(gene=human.pos.MCL[,1],
                             keyType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 1)

human.neg.MCL.go <- enrichGO(gene=human.neg.MCL[,1],
                             keyType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 1)

chimp.pos.MCL.go <- enrichGO(gene=chimp.pos.MCL[,1],
                             keyType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 1)

chimp.neg.MCL.go <- enrichGO(gene=chimp.neg.MCL[,1],
                             keyType = "SYMBOL", 
                             OrgDb = org.Hs.eg.db, 
                             ont = "BP", 
                             pAdjustMethod = "BH",
                             pvalueCutoff = 0.05, 
                             qvalueCutoff = 1)

dotplot(human.pos.MCL.go, showCategory=15) + ggtitle("human.pos-specific cluster")
dotplot(human.neg.MCL.go, showCategory=15) + ggtitle("human.neg-specific cluster")
dotplot(chimp.pos.MCL.go, showCategory=15) + ggtitle("chimp.pos-specific cluster")
dotplot(chimp.neg.MCL.go, showCategory=15) + ggtitle("chimp.neg-specific cluster")

#dotplot
MCL1 <- dotplot(human.pos.MCL.go, showCategory=15) + ggtitle("human.pos-specific cluster")
MCL1 <- MCL1 + scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none') 
MCL2 <- dotplot(human.neg.MCL.go, showCategory=15) + ggtitle("human.neg-specific cluster")
MCL2 <- MCL2+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
MCL3 <- dotplot(chimp.pos.MCL.go, showCategory=15) + ggtitle("chimp.pos-specific cluster")
MCL3 <- MCL3+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
MCL4 <- dotplot(chimp.neg.MCL.go, showCategory=15) + ggtitle("chimp.neg-specific cluster")
MCL4 <- MCL4+ scale_y_discrete(labels= label_wrap(30))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')


cowplot::plot_grid(MCL1, MCL2, MCL3, MCL4, ncol=2, labels=LETTERS[1:4])
#plot size 900x900

#emapplot
MCL_merg <- rbind(human.pos.MCL, human.neg.MCL, chimp.pos.MCL, chimp.neg.MCL)
write.table(MCL_merg, file = "correlatedgenes_species-specific_MCL.txt", row.names = F, quote = F)

MCL_merg_res <- compareCluster(ID~group, data=MCL_merg, fun="enrichGO", OrgDb = org.Hs.eg.db,keyType = "SYMBOL", pvalueCutoff=0.05)
MCL_merg_res <- pairwise_termsim(MCL_merg_res)
emapplot(MCL_merg_res) 
