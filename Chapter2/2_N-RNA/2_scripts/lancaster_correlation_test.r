#############################################################################
#correlation test with merged read count table 
#ZEB2 expression vs all other gene expression during 7 time points

#add gene symbols to ensembl id
#Vladi would like to have both information
#############################################################################
#library calling
library(psych)
library(dplyr)

#table with normalized read counts
cts_norm

#Add gene symbols
ensembl <- useMart("ensembl")
listDatasets(ensembl)
ensembl.human <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
attributes <- listAttributes(ensembl.human)
filters = listFilters(ensembl.human)

symbol = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
               filters = 'ensembl_gene_id', 
               values = rownames(cts_norm), 
               mart = ensembl.human)
nrow(symbol) #41014 genes
nrow(cts_norm) #41122 genes

#sort the gene id in symbol data according to the rownames of data
symbol.sort <- symbol[match(rownames(cts_norm), symbol$ensembl_gene_id),]
#add a new column of gene symbols sorted by rownames of data
cts_norm_symbol <- cbind(symbol.sort$external_gene_name, cts_norm)
colnames(cts_norm_symbol)[1] <- "gene symbol"

#run correlation test
tmp_cor <- data.frame(matrix(, nrow=nrow(cts_norm), ncol=1))
rownames(tmp_cor) <- rownames(cts_norm)
colnames(tmp_cor) <- "cor"

#pearson correlation test
corr_h <-data.frame(matrix(, nrow=nrow(cts_norm), ncol=7))
rownames(corr_h) <- rownames(cts_norm)
colnames(corr_h) <- h

for (i in 1:length(h)){
  tmp_stage <- as.matrix(cts_norm[,which(stringr::str_detect(Samples_h, paste(h[i], "_", sep = "")))])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554.21"),]), as.numeric(tmp_stage[j,]),method ="pearson", adjust="BH")$p.adj
    tmp_cor[is.na(tmp_cor)] <- 1
    if (tmp_cor[j,] < 0.05){
      corr_h[j, h[i]] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554.21"),]), as.numeric(tmp_stage[j,]),method ="pearson")$r
    } 
  }
}

#remove ZEB2 row
corr_h <- corr_h[-which(rownames(corr_h) =="ENSG00000169554.21"),]

#spearman correlation test

#corr_h_s <-data.frame(matrix(, nrow=nrow(cts_norm_symbol), ncol=8))
#rownames(corr_h_s) <- rownames(cts_norm_symbol)
#corr_h_s[,1] <- cts_norm_symbol[,1]
#colnames(corr_h_s) <- c("gene symbol", h)

#for (i in 1:length(h)){
#  tmp_stage <- as.matrix(cts_norm_symbol[,which(stringr::str_detect(Samples_h, paste(h[i], "_", sep = "")))+1])
#  for (j in 1:nrow(tmp_stage)){
#    tmp_cor[j,] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554"),]), as.numeric(tmp_stage[j,]),method ="spearman", adjust="BH")$p.adj
#    tmp_cor[is.na(tmp_cor)] <- 1
#    if (tmp_cor[j,] < 0.05){
#      corr_h_s[j, h[i]] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554"),]), as.numeric(tmp_stage[j,]),method ="spearman")$r
#    } 
#  }
#}


#############################################################################
#test with p adj
#############################################################################
summary(attitude)
head(attitude)
class(attitude)

ct <- corr.test(attitude)  #find the correlations and give the probabilities
ct #show the results
cts <- corr.test(attitude[1],attitude[4:6]) #reports all values corrected for multiple tests
cts$p.adj
cts_adj <- corr.test(attitude[1:3],attitude[4:6],adjust="BH")
cts_adj$p.adj
#corr.test(sat.act[1:3],sat.act[4:6],adjust="none")  #don't adjust the probabilities

#take correlations and show the probabilities as well as the confidence intervals
print(corr.p(cts$r,n=30),short=FALSE)  
x <- corr.p(cts$r,n=30)
x$p
y <- corr.p(cts_adj$r,n=30)
y$p
library(rstatix)
adjust_pvalue(cts, method = "BH")


tmp_stage_x <- lapply(tmp_stage,as.numeric)

test_corr <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554"),]), as.numeric(tmp_stage[56,]),method ="pearson", adjust="BH")
test_corr_v2 <-corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554"),]), as.data.frame(tmp_stage[-5523,]),method ="pearson", adjust="BH")
test_corr$n
test_corr$p.adj
print(corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554"),]), as.numeric(tmp_stage[56,]),method ="pearson", adjust="BH"))

corr.p(tmp_cor$cor, n=nrow(tmp_cor), adjust="BH", alpha =.05)
#############################################################################
tmp_cor <- data.frame(matrix(, nrow=nrow(cts_norm), ncol=1))
rownames(tmp_cor) <- rownames(cts_norm)
colnames(tmp_cor) <- "cor"

#pearson correlation test
corr_g <-data.frame(matrix(, nrow=nrow(cts_norm), ncol=7))
rownames(corr_g) <- rownames(cts_norm)
colnames(corr_g) <- g

for (i in 1:length(g)){
  tmp_stage <- as.matrix(cts_norm[,which(stringr::str_detect(Samples_g, paste(g[i], "_", sep = "")))+21])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554.21"),]), as.numeric(tmp_stage[j,]),method ="pearson", adjust="BH")$p.adj
    tmp_cor[is.na(tmp_cor)] <- 1
    if (tmp_cor[j,] < 0.05){
      corr_g[j, g[i]] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554.21"),]), as.numeric(tmp_stage[j,]),method ="pearson")$r
    } 
  }
}

#remove ZEB2 row
corr_g <- corr_g[-which(rownames(corr_g) =="ENSG00000169554.21"),]

#spearman correlation test

#corr_g_s <-data.frame(matrix(, nrow=nrow(cts_norm_symbol), ncol=8))
#rownames(corr_g_s) <- rownames(cts_norm_symbol)
#corr_g_s[,1] <- cts_norm_symbol[,1]
#colnames(corr_g_s) <- c("gene symbol", g)

#for (i in 1:length(g)){
#  tmp_stage <- as.matrix(cts_norm_symbol[,which(stringr::str_detect(Samples_g, paste(g[i], "_", sep = "")))+22])
#  for (j in 1:nrow(tmp_stage)){
#    tmp_cor[j,] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554"),]), as.numeric(tmp_stage[j,]),method ="spearman", adjust="BH")$p.adj
#    tmp_cor[is.na(tmp_cor)] <- 1
#    if (tmp_cor[j,] < 0.05){
#      corr_g_s[j, g[i]] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554"),]), as.numeric(tmp_stage[j,]),method ="spearman")$r
#    } 
#  }
#}

#correlated genes with ZEB2 at day 0 + after re-version
#the final number of correlated genes of humans stays the same, but of gorilla changes back to the first try
length(na.omit(corr_h$H9_d0)) #1649 #1650 #the final one w/o ZEB2 1649
length(na.omit(corr_g$G1_d0)) #2709 #4247 #the final one w/o ZEB2 2709

#correlated genes with ZEB2 at day 2
length(na.omit(corr_h$H9_d2)) #1610 #1611 #the final one w/o ZEB2 1610
length(na.omit(corr_g$G1_d2)) #1769 #3285 #the final one w/o ZEB2 1768

#correlated genes with ZEB2 at day 3
length(na.omit(corr_h$H9_d3)) #2579 #2581 #the final one w/o ZEB2 2580
length(na.omit(corr_g$G1_d3)) #1251 #3766 #the final one w/o ZEB2 1250

#correlated genes with ZEB2 at day 5
length(na.omit(corr_h$H9_d5)) #1572 #1572 #the final one w/o ZEB2 1571
length(na.omit(corr_g$G1_d5)) #1074 #2575 #the final one w/o ZEB2 1073

#correlated genes with ZEB2 at day 10
length(na.omit(corr_h$H9_d10)) #1323 #1323 #the final one w/o ZEB2 1322
length(na.omit(corr_g$G1_d10)) #1710 #2970 #the final one w/o ZEB2 1709

#correlated genes with ZEB2 at day 15
length(na.omit(corr_h$H9_d15)) #2064 #2064 #the final one w/o ZEB2 2063
length(na.omit(corr_g$G1_d15)) #1843 #3796 #the final one w/o ZEB2 1842

#correlated genes with ZEB2 at day 25
length(na.omit(corr_h$H9_d25)) #1033 #1033 #the final one w/o ZEB2 1032
length(na.omit(corr_g$G1_d25)) #1276 #2257 #the final one w/o ZEB2 1275

saveRDS(corr_h, file = here("output/preliminary/corr_h_FINAL.rds"))
saveRDS(corr_g, file = here("output/preliminary/corr_g_FINAL.rds"))

path_out <- "output/correlation_GTF/"

for (i in 1:ncol(corr_h)) {
  corr_h_tmp <- data.frame("ensembl ID" = rownames(corr_h)[(!is.na(corr_h[,i]))], 
                           "gene symbol" =ID_symbol[which(ID_symbol$Geneid %in% rownames(corr_h)[which(!is.na(corr_h[,i]))]),2])
  write.table(corr_h_tmp, file = paste(path_out,'correlatedGeneswithZEB2_FINAL_',colnames(corr_h)[i],'_n',length(na.omit(corr_h[,i])),'.tsv', sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
}


for (i in 1:ncol(corr_g)) {
  corr_g_tmp <- data.frame("ensembl ID" = rownames(corr_g)[(!is.na(corr_g[,i]))], 
                           "gene symbol" =ID_symbol[which(ID_symbol$Geneid %in% rownames(corr_g)[which(!is.na(corr_g[,i]))]),2])
  write.table(corr_g_tmp, file = paste(path_out,'correlatedGeneswithZEB2_FINAL_',colnames(corr_g)[i],'_n',length(na.omit(corr_g[,i])),'.txt', sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
}
#ZEB2 exists in all 14 lists above

#############################################################################
#overlapping genes at one time point (human vs gorilla)
#############################################################################

length(rownames(cts_norm_symbol)[(which(rownames(corr_h)[(!is.na(corr_h[,2]))] %in% rownames(corr_g)[(!is.na(corr_g[,2]))]))]) #113
length(cts_norm_symbol[which(rownames(cts_norm_symbol) %in% rownames(corr_h)[(!is.na(corr_h[,2]))][which(rownames(corr_h)[(!is.na(corr_h[,2]))] %in% rownames(corr_g)[(!is.na(corr_g[,2]))])])]) # also 113

length(rownames(cts_norm_symbol)[(which(rownames(corr_h)[(!is.na(corr_h[,3]))] %in% rownames(corr_g)[(!is.na(corr_g[,3]))]))]) #95

length(rownames(cts_norm_symbol)[(which(rownames(corr_h)[(!is.na(corr_h[,4]))] %in% rownames(corr_g)[(!is.na(corr_g[,4]))]))]) #66

length(rownames(cts_norm_symbol)[(which(rownames(corr_h)[(!is.na(corr_h[,5]))] %in% rownames(corr_g)[(!is.na(corr_g[,5]))]))]) #71

length(rownames(cts_norm_symbol)[(which(rownames(corr_h)[(!is.na(corr_h[,6]))] %in% rownames(corr_g)[(!is.na(corr_g[,6]))]))]) #63

length(rownames(cts_norm_symbol)[(which(rownames(corr_h)[(!is.na(corr_h[,7]))] %in% rownames(corr_g)[(!is.na(corr_g[,7]))]))]) #111

length(rownames(cts_norm_symbol)[(which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,8]))]))]) #52

for (i in 2:ncol(corr_h)) {
  corr_overlap_tp_tmp <- as.data.frame(cts_norm_symbol[which(rownames(cts_norm_symbol) %in% rownames(corr_h)[(!is.na(corr_h[,i]))][which(rownames(corr_h)[(!is.na(corr_h[,i]))] %in% rownames(corr_g)[(!is.na(corr_g[,i]))])]),1])
  write.table(corr_overlap_tp_tmp, 
              file = paste(path_out,
                           'overlapping_corrgenes_HvsG',
                           '_d',
                           x[i-1],
                           '_n',
                           length(rownames(cts_norm_symbol)[(which(rownames(corr_h)[(!is.na(corr_h[,i]))] %in% rownames(corr_g)[(!is.na(corr_g[,i]))]))]),
                           '.txt', 
                           sep = ""), 
              row.names = T, col.names = F)
}

#ZEB2 exists in all 7 lists above

#overlapps with pos/neg
sum(rownames(corr_h)[which(sign(corr_h$H9_d0)==1)] %in% rownames(corr_g)[which(sign(corr_g$G1_d0)==1)]) #+30

sum(rownames(corr_h)[which(sign(corr_h$H9_d0)==-1)] %in% rownames(corr_g)[which(sign(corr_g$G1_d0)==-1)]) #-22

sum(rownames(corr_h)[which(sign(corr_h$H9_d0)==1)] %in% rownames(corr_g)[which(sign(corr_g$G1_d0)==-1)]) #+/-34

sum(rownames(corr_h)[which(sign(corr_h$H9_d0)==-1)] %in% rownames(corr_g)[which(sign(corr_g$G1_d0)==1)]) #-/+27

cts_norm_symbol[which(rownames(cts_norm_symbol) %in% rownames(corr_h)[which(sign(corr_h[,2])==1)][which(rownames(corr_h)[which(sign(corr_h[,2])==1)] %in% rownames(corr_g)[which(sign(corr_g[,2])==1)])]),1]

for (i in 2:ncol(corr_h)) {
  corr_overlap_tp_pos_tmp <- as.data.frame(cts_norm_symbol[which(rownames(cts_norm_symbol) %in% rownames(corr_h)[which(sign(corr_h[,i])==1)][which(rownames(corr_h)[which(sign(corr_h[,i])==1)] %in% rownames(corr_g)[which(sign(corr_g[,i])==1)])]),1])
  colnames(corr_overlap_tp_pos_tmp)[1] <- "gene_symbol"
  corr_overlap_tp_pos_tmp$corr <- "(+)"
  corr_overlap_tp_neg_tmp <- as.data.frame(cts_norm_symbol[which(rownames(cts_norm_symbol) %in% rownames(corr_h)[which(sign(corr_h[,i])==-1)][which(rownames(corr_h)[which(sign(corr_h[,i])==-1)] %in% rownames(corr_g)[which(sign(corr_g[,i])==-1)])]),1])
  colnames(corr_overlap_tp_neg_tmp)[1] <- "gene_symbol"
  corr_overlap_tp_neg_tmp$corr <- "(-)"
  corr_overlap_tp_posH_tmp <- as.data.frame(cts_norm_symbol[which(rownames(cts_norm_symbol) %in% rownames(corr_h)[which(sign(corr_h[,i])==1)][which(rownames(corr_h)[which(sign(corr_h[,i])==1)] %in% rownames(corr_g)[which(sign(corr_g[,i])==-1)])]),1])
  colnames(corr_overlap_tp_posH_tmp)[1] <- "gene_symbol"
  corr_overlap_tp_posH_tmp$corr <- "h(+)g(-)"
  corr_overlap_tp_negH_tmp <- as.data.frame(cts_norm_symbol[which(rownames(cts_norm_symbol) %in% rownames(corr_h)[which(sign(corr_h[,i])==-1)][which(rownames(corr_h)[which(sign(corr_h[,i])==-1)] %in% rownames(corr_g)[which(sign(corr_g[,i])==1)])]),1])
  colnames(corr_overlap_tp_negH_tmp)[1] <- "gene_symbol"
  corr_overlap_tp_negH_tmp$corr <- "h(-)g(+)"
  corr_overlap_tp_p_n_tmp <- rbind(corr_overlap_tp_pos_tmp, corr_overlap_tp_neg_tmp, corr_overlap_tp_posH_tmp, corr_overlap_tp_negH_tmp) 
  write.table(corr_overlap_tp_p_n_tmp, 
              file = paste(path_out,
                           'overlapping_corrgenes_HvsG_pos_neg',
                           '_d',
                           x[i-1],
                           '.txt', 
                           sep = ""), 
              row.names = T, col.names = F)
}


#finding overlapping genes disregarding time points, e.g. gorilla day 3 vs human day 10
#overlaps between all pairs of time points and represent this as a matrix

as.data.frame(cts_norm_symbol[which(rownames(cts_norm_symbol) %in% rownames(corr_h)[(!is.na(corr_h[,1]))][which(rownames(corr_h)[(!is.na(corr_h[,i]))] %in% rownames(corr_g)[(!is.na(corr_g[,i]))])]),1])

#test with igraph
library(igraph)

nr_overlaps <- c()
links_col<- c()
link_tmp <- NULL

for (i in 2:ncol(corr_h)){
  nr_overlaps[length(nr_overlaps) + 1] <- c(length(rownames(corr_h)[(!is.na(corr_h[,8]))][which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,i]))])]))
  links_col <- c(links_col, rownames(corr_h)[(!is.na(corr_h[,8]))][which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,i]))])])
  link_tmp <- data.frame(matrix(, nrow=7, ncol=length(unique(links_col))))
  rownames(link_tmp) <- paste("day", x, sep = "")
  colnames(link_tmp) <- unique(links_col)
}

overlap_tbl  <- data.frame(matrix(, nrow=7, ncol=7))
overlap_tbl$X4 <- nr_overlaps
colnames(overlap_tbl) <- paste("H_day", x, sep = "")
rownames(overlap_tbl) <- paste("G_day", x, sep = "")

path_out <- here("output/correlation/")
write.table(overlap_tbl, file = paste(path_out,'overlapping_genes_HvsG_overallcomparison','.txt', sep = ""), row.names = T, col.names = T)


link_tmp[1,] <- ifelse(colnames(link_tmp[1,]) %in% rownames(corr_h)[(!is.na(corr_h[,8]))][which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,2]))])],1, 0)
link_tmp[2,] <- ifelse(colnames(link_tmp[2,]) %in% rownames(corr_h)[(!is.na(corr_h[,8]))][which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,3]))])],1, 0)
link_tmp[3,] <- ifelse(colnames(link_tmp[3,]) %in% rownames(corr_h)[(!is.na(corr_h[,8]))][which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,4]))])],1, 0)
link_tmp[4,] <- ifelse(colnames(link_tmp[4,]) %in% rownames(corr_h)[(!is.na(corr_h[,8]))][which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,5]))])],1, 0)
link_tmp[5,] <- ifelse(colnames(link_tmp[5,]) %in% rownames(corr_h)[(!is.na(corr_h[,8]))][which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,6]))])],1, 0)
link_tmp[6,] <- ifelse(colnames(link_tmp[6,]) %in% rownames(corr_h)[(!is.na(corr_h[,8]))][which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,7]))])],1, 0)
link_tmp[7,] <- ifelse(colnames(link_tmp[7,]) %in% rownames(corr_h)[(!is.na(corr_h[,8]))][which(rownames(corr_h)[(!is.na(corr_h[,8]))] %in% rownames(corr_g)[(!is.na(corr_g[,8]))])],1, 0)

dim(link_tmp)

nodes <- NULL
nodes <- data.frame(matrix(, nrow= nrow(link_tmp)+ncol(link_tmp), ncol=2))
colnames(nodes) <- c("id", "days")
nodes[,1] <- c(rownames(link_tmp), colnames(link_tmp))

nodes[1,2] <- "day0"
nodes[2,2] <- "day2"
nodes[3,2] <- "day3"
nodes[4,2] <- "day5"
nodes[5,2] <- "day10"
nodes[6,2] <- "day15"
nodes[7,2] <- "day25"

net_H0 <- graph_from_incidence_matrix(link_tmp)


# Generate colors based on days
V(net_H0)$color <- c("steel blue", "orange")[V(net_H0)$type+1]
colrs <- c("#fa938e", "#d0af47", "#7bc445", "#52ccab", "#50c4ef", "#b9a2ff", "#fb82df")
V(net_H0)$color[V(net_H0)$type==F] <- colrs[V(net_H0)$type==F] 
V(net_H0)$color[V(net_H0)$name =="ENSG00000169554"] <- "red" #mark ZEB2 red

# genes uniquely overlapping at a certain time point = same color as the days
for (i in 1:nrow(link_tmp)){
  unique <- setdiff(colnames(link_tmp)[which(link_tmp[i,] == "1")], colnames(link_tmp[colSums((link_tmp)) > 1]))
  V(net_H0)$color[which(V(net_H0)$name %in% unique)] <- colrs[i]
}


# days will have name labels, genes will not:
V(net_H0)$label <- ""
V(net_H0)$label[V(net_H0)$type==F] <- nodes$days[V(net_H0)$type==F] 
V(net_H0)$label.cex=6
V(net_H0)$label.font=2

order <- c("day0", links_col[1:nr_overlaps[1]])
order <- c(order, "day2", links_col[nr_overlaps[1]:sum(nr_overlaps[1:2])])
order <- c(order, "day3", links_col[sum(nr_overlaps[1:2]):sum(nr_overlaps[1:3])])
order <- c(order, "day5", links_col[sum(nr_overlaps[1:3]):sum(nr_overlaps[1:4])])
order <- c(order, "day10", links_col[sum(nr_overlaps[1:4]):sum(nr_overlaps[1:5])])
order <- c(order, "day15", links_col[sum(nr_overlaps[1:5]):sum(nr_overlaps[1:6])])
order <- c(order, "day25", links_col[sum(nr_overlaps[1:6]):sum(nr_overlaps[1:7])])

coords <- layout_in_circle(net_H0, order = order)

plot(net_H0, layout=coords, vertex.label.color="black", vertex.size=(2-V(net_H0)$type)*4)

#networks are saved as 900x550
#circular networks are saved as 800x800

#############################################################################
#overlapping genes at all time points (human vs gorilla)
#############################################################################
h_day0 <- corr_h %>% dplyr::select(2) %>% na.omit %>% rownames()
h_day2 <- corr_h %>% dplyr::select(3) %>% na.omit %>% rownames()
h_day3 <- corr_h %>% dplyr::select(4) %>% na.omit %>% rownames()
h_day5 <- corr_h %>% dplyr::select(5) %>% na.omit %>% rownames()
h_day10 <- corr_h %>% dplyr::select(6) %>% na.omit %>% rownames()
h_day15 <- corr_h %>% dplyr::select(7) %>% na.omit %>% rownames()
h_day25 <- corr_h %>% dplyr::select(8) %>% na.omit %>% rownames()

g_day0 <- corr_g %>% dplyr::select(2) %>% na.omit %>% rownames()
g_day2 <- corr_g %>% dplyr::select(3) %>% na.omit %>% rownames()
g_day3 <- corr_g %>% dplyr::select(4) %>% na.omit %>% rownames()
g_day5 <- corr_g %>% dplyr::select(5) %>% na.omit %>% rownames()
g_day10 <- corr_g %>% dplyr::select(6) %>% na.omit %>% rownames()
g_day15 <- corr_g %>% dplyr::select(7) %>% na.omit %>% rownames()
g_day25 <- corr_g %>% dplyr::select(8) %>% na.omit %>% rownames()

#############################################################################
#correlation test before and after ZEB2 in each species
#############################################################################
#human
corr_h_peak <-data.frame(matrix(, nrow=nrow(cts_norm), ncol=2))
rownames(corr_h_peak) <- rownames(cts_norm)
colnames(corr_h_peak) <- c("pre", "post")

tmp_stage_e <- as.matrix(cts_norm[,which(stringr::str_detect(Samples_h, "d0_|d2_|d3_|d5_|d10_"))])
tmp_stage_l <- as.matrix(cts_norm[,which(stringr::str_detect(Samples_h, "d10_|d15_|d25_"))])

for (i in 1:nrow(tmp_stage_e)){
  tmp_cor[i,] <- corr.test(as.numeric(tmp_stage_e[which(rownames(tmp_stage_e) =="ENSG00000169554.21"),]), as.numeric(tmp_stage_e[i,]),method ="pearson", adjust="BH")$p.adj
  tmp_cor[is.na(tmp_cor)] <- 1
    if (tmp_cor[i,] < 0.05){
      corr_h_peak[i, 1] <- corr.test(as.numeric(tmp_stage_e[which(rownames(tmp_stage_e) =="ENSG00000169554.21"),]), as.numeric(tmp_stage_e[i,]),method ="pearson")$r
    }
}


for (j in 1:nrow(tmp_stage_l)){
  tmp_cor[j,] <- corr.test(as.numeric(tmp_stage_l[which(rownames(tmp_stage_l) =="ENSG00000169554.21"),]), as.numeric(tmp_stage_l[j,]),method ="pearson", adjust="BH")$p.adj
  tmp_cor[is.na(tmp_cor)] <- 1
  if (tmp_cor[j,] < 0.05){
    corr_h_peak[j, 2] <- corr.test(as.numeric(tmp_stage_l[which(rownames(tmp_stage_l) =="ENSG00000169554.21"),]), as.numeric(tmp_stage_l[j,]),method ="pearson")$r
  }
}

nrow(corr_h_peak) #41204
nrow(na.omit(corr_h_peak)) #3713


#gorilla
corr_g_peak <-data.frame(matrix(, nrow=nrow(cts_norm), ncol=2))
rownames(corr_g_peak) <- rownames(cts_norm)
colnames(corr_g_peak) <- c("pre", "post")

tmp_stage_g_e <- as.matrix(cts_norm[,which(stringr::str_detect(Samples_g, "d0_|d2_|d3_|d5_"))+1])
tmp_stage_g_l <- as.matrix(cts_norm[,which(stringr::str_detect(Samples_g, "d5_|d10_|d15_|d25_"))+1])


for (i in 1:nrow(tmp_stage_g_e)){
  tmp_cor[i,] <- corr.test(as.numeric(tmp_stage_g_e[which(rownames(tmp_stage_g_e) =="ENSG00000169554.21"),]), as.numeric(tmp_stage_g_e[i,]),method ="pearson", adjust="BH")$p.adj
  tmp_cor[is.na(tmp_cor)] <- 1
  if (tmp_cor[i,] < 0.05){
    corr_g_peak[i, 1] <- corr.test(as.numeric(tmp_stage_g_e[which(rownames(tmp_stage_g_e) =="ENSG00000169554.21"),]), as.numeric(tmp_stage_g_e[i,]),method ="pearson")$r
  }
}

  
  
for (j in 1:nrow(tmp_stage_g_l)){
  tmp_cor[j,] <- corr.test(as.numeric(tmp_stage_g_l[which(rownames(tmp_stage_g_l) =="ENSG00000169554.21"),]), as.numeric(tmp_stage_g_l[j,]),method ="pearson", adjust="BH")$p.adj
  tmp_cor[is.na(tmp_cor)] <- 1
  if (tmp_cor[j,] < 0.05){
    corr_g_peak[j, 2] <- corr.test(as.numeric(tmp_stage_g_l[which(rownames(tmp_stage_g_l) =="ENSG00000169554.21"),]), as.numeric(tmp_stage_g_l[j,]),method ="pearson")$r
  }
}

#############################################################################
#correlation test throughout the time range
#############################################################################
tmp_cor <- data.frame(matrix(, nrow=nrow(cts_norm), ncol=1))
rownames(tmp_cor) <- rownames(cts_norm)
colnames(tmp_cor) <- "cor"

corr_h_ts <-data.frame(matrix(, nrow=nrow(cts_norm), ncol=1))
rownames(corr_h_ts) <- rownames(cts_norm)
colnames(corr_h_ts) <- c("coef")

for (i in 1:length(h)){
  tmp_stage <- as.matrix(cts_norm[,1:21])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554.21"),]), as.numeric(tmp_stage[j,]),method ="pearson", adjust="BH")$p.adj
    tmp_cor[is.na(tmp_cor)] <- 1
    if (tmp_cor[j,] < 0.05){
      corr_h_ts[j, 1] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554.21"),]), as.numeric(tmp_stage[j,]),method ="pearson")$r
    } 
  }
}

nrow(corr_h_ts) #41204
corr_h_ts <- na.omit(corr_h_ts) 
nrow(corr_h_ts) #10232 this contains ZEB2
saveRDS(corr_h_ts, file = here("output/preliminary/corr_h_ts_FINAL.rds"))

corr_h_ts <- readRDS(file = "output/preliminary/corr_h_ts_FINAL.rds")

corr_g_ts <-data.frame(matrix(, nrow=nrow(cts_norm), ncol=1))
rownames(corr_g_ts) <- rownames(cts_norm)
colnames(corr_g_ts) <- c("coef")

for (i in 1:length(g)){
  tmp_stage <- as.matrix(cts_norm[,22:42])
  for (j in 1:nrow(tmp_stage)){
    tmp_cor[j,] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554.21"),]), as.numeric(tmp_stage[j,]),method ="pearson", adjust="BH")$p.adj
    tmp_cor[is.na(tmp_cor)] <- 1
    if (tmp_cor[j,] < 0.05){
      corr_g_ts[j, 1] <- corr.test(as.numeric(tmp_stage[which(rownames(tmp_stage) =="ENSG00000169554.21"),]), as.numeric(tmp_stage[j,]),method ="pearson")$r
    } 
  }
}

nrow(corr_g_ts) #41204
corr_g_ts <- na.omit(corr_g_ts) 
nrow(corr_g_ts) #6293  this contains ZEB2
saveRDS(corr_g_ts, file = here("output/preliminary/corr_g_ts_FINAL.rds"))

corr_g_ts <- readRDS(file = "output/preliminary/corr_g_ts_FINAL.rds")

#############################################################################
#cutoff of |r| < 0.5
#############################################################################
#Human result table after cutoff of |r| < 0.5
nrow(corr_h_ts) #10232
corr_h_ts_0.5 <- corr_h_ts
corr_h_ts_0.5[abs(corr_h_ts_0.5) < 0.5] <- NA
corr_h_ts_0.5 <- na.omit(corr_h_ts_0.5)
nrow(corr_h_ts_0.5) #7508

#Gorilla result table after cutoff of |r| < 0.5
nrow(corr_g_ts) #6293
corr_g_ts_0.5 <- corr_g_ts
corr_g_ts_0.5[abs(corr_g_ts_0.5) < 0.5] <- NA
corr_g_ts_0.5 <- na.omit(corr_g_ts_0.5)
nrow(corr_g_ts_0.5) #4025

#distribution of r values
corr_h_ts$cond <- "BEFORE"
corr_h_ts_0.5$cond <- "AFTER"
corr_h_hist <- rbind(corr_h_ts, corr_h_ts_0.5)
human_hist <- ggplot(corr_h_hist, aes(coef, fill = cond)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') +
  ggtitle("Distribution of r-values (human)")

corr_g_ts$cond <- "BEFORE"
corr_g_ts_0.5$cond <- "AFTER"
corr_g_hist <- rbind(corr_g_ts, corr_g_ts_0.5)
gorilla_hist <- ggplot(corr_g_hist, aes(coef, fill = cond)) + 
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity') +
  ggtitle("Distribution of r-values (gorilla)")

ggarrange(human_hist,
          gorilla_hist, 
          ncol = 1, nrow = 2)

#data output
corr_h_ts_0.5$genename <- ID_symbol %>% filter(Geneid %in% rownames(corr_h_ts_0.5)) %>% pull(-1)
corr_g_ts_0.5$genename <- ID_symbol %>% filter(Geneid %in% rownames(corr_g_ts_0.5)) %>% pull(-1)

path_out <- "../../Volumes/Extreme SSD/PhD/Drylab/Lancaster/output/correlation_GTF/"
 
write.table(corr_h_ts_0.5, 
            file = paste(path_out,'correlatedGeneswithZEB2_FINAL_cutoff0.5_alltime_human_', nrow(corr_h_ts_0.5), '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

write.table(corr_g_ts_0.5, 
            file = paste(path_out,'correlatedGeneswithZEB2_FINAL_cutoff0.5_alltime_gorilla_', nrow(corr_g_ts_0.5), '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

#both gene lists contain ZEB2
