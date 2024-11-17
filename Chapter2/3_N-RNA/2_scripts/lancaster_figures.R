#############################################################################
#Plotting with merged read count table - gene expression, e.g. ZEB2/PCA
#############################################################################
#as same way how I did in B-cell RNA-seq analysis
cts_norm %>% as.data.frame() %>% filter(rownames(cts_norm) == "ENSG00000169554.21") %>%
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>%
  mutate(time = rep(tp, 2)) %>%
  mutate(species = c(rep("Human", 21), rep("Gorilla", 21))) %>%
  mutate(time = factor(time, levels = stri_paste("day",x))) %>%
  mutate(species = factor(species, levels = c("Human", "Gorilla"))) %>%
  ggplot(aes(x = time, y = count, color = species, group = time)) +
  geom_boxplot() + stat_summary(fun.y = mean, geom = "line") +
  scale_color_manual(values=c("#8dd35fff", "#73C2FB")) +
  scale_y_log10() +
  facet_wrap( ~ species) + theme_bw() + labs(title="Counts of ZEB2 in each species", x="") +
  ylab("normalized count (log)")

#no facet
cts_norm %>% as.data.frame() %>% filter(rownames(cts_norm) == "ENSG00000169554.21") %>%
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>%
  mutate(time = rep(tp, 2)) %>%
  mutate(species = c(rep("Human", 21), rep("Gorilla", 21))) %>%
  mutate(time = factor(time, levels = stri_paste("day",x))) %>%
  mutate(species = factor(species, levels = c("Human", "Gorilla"))) %>%
  ggplot(aes(x = time, y = count, fill = species)) +
  geom_boxplot() + stat_summary(fun.y = mean, geom = "line") +
  scale_fill_manual(values=c("#8dd35fff", "#73C2FB")) +
  scale_y_log10() + geom_dotplot(binaxis='y', stackdir='center',
                                 position=position_dodge(0.75), dotsize=0.5) +
  theme_bw() + labs(title="Counts of ZEB2 in each species", x="") +
  ylab("normalized count (log)")

#plot the counts for the groups over time for the gene with the smallest adjusted p value, 
#testing for condition-dependent time profile and accounting for differences at time 0

#Plot the counts for the groups over time for ZEB2
plot.tp.sp_ZEB2 <- plotCounts(dds, which(rownames(res) =="ENSG00000169554.21"), 
                              intgroup = c("tp","sp"), returnData = TRUE)

ggplot(plot.tp.sp_ZEB2,
       aes(x = tp, y = count, color = sp, group = tp)) + 
  geom_point() + stat_summary(fun.y=mean, geom="line") +
  scale_y_log10() +
  labs(title="Counts for species over time for ZEB2", x="") +
  guides(color=guide_legend("species")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(plot.title = element_text(face="bold", size = 16)) +
  geom_point(size = 2) +
  scale_color_manual(values=c("#8dd35fff", "#73C2FB")) +
  facet_wrap( ~ sp) +
  stat_summary(fun = mean, geom = "line") 

#Size of ggplot above: 550x530

#PCA
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup=c("sp")) + #PCA with species
  guides(color=guide_legend("species")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(plot.title = element_text(face="bold", size = 16)) +
  geom_point(size = 2) 

plotPCA(vsd, intgroup=c("tp")) + #PCA with timepoints
  guides(color=guide_legend("timpoints")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(plot.title = element_text(face="bold", size = 16)) +
  geom_point(size = 2) 

pcaData <- plotPCA(vsd, intgroup=c("sp", "tp"), returnData = TRUE) #PCA with species and tp together
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color = tp, shape = sp)) +
  geom_point(size = 3) +
  theme_bw() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) 

#Size of ggplot above: 550x530

# Modified plotPCA from DESeq2 package to show pc3
library(genefilter)
library(ggplot2)
library(ggrepel)

plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(PC2 = pca$x[, 2], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(object)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[2:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC2", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC2: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}


plotPCA.san(vsd, intgroup=c("sp")) + #PCA with species
  guides(color=guide_legend("species")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(plot.title = element_text(face="bold", size = 16)) +
  geom_point(size = 2) 

#############################################################################
#Plotting the number of correlated genes with ZEB2
#############################################################################
#a matrix with h time & g time showing the number overlapping genes (pos & neg at each time point)
nr_corr <- data.frame(matrix(, nrow=length(h)*4, ncol=4))
colnames(nr_corr) <- c("timepoint", "correlation", "species","n")
nr_corr[,1] <- paste("day", x, sep = "")
nr_corr[,2] <- c(rep("pos", 7), rep("neg",7))
nr_corr[,3] <- c(rep("human", length(h)*2), rep("gorilla", length(h)*2))

for (i in 1:nrow(nr_corr)){
  if (i < 15) {
    if (i < 8){
      nr_corr[i,4] <- as.numeric(length(rownames(corr_h)[which(sign(corr_h[,i])==1)]))
    } else {
      nr_corr[i,4] <- as.numeric(length(rownames(corr_h)[which(sign(corr_h[,i-7])==-1)]))
    }
  } else {
    if (i < 22) {
      nr_corr[i,4] <- as.numeric(length(rownames(corr_g)[which(sign(corr_g[,i-14])==1)]))
    } else {
      nr_corr[i,4] <- as.numeric(length(rownames(corr_g)[which(sign(corr_g[,i-21])==-1)]))
    }
  }
}

#change the nr of corr into numeric
nr_corr$n <- as.numeric(nr_corr$n)
#set the order 
nr_corr$timepoint <- factor(nr_corr$timepoint,levels = c("day0", "day2", "day3", "day5", "day10", "day15", "day25"))
nr_corr$species <- factor(nr_corr$species, levels = c("human", "gorilla"))
nr_corr$correlation <- factor(nr_corr$correlation, levels = c("pos", "neg"))


library(plyr)
# Sort by timpoint and species
nr_corr_sort <- arrange(nr_corr, timepoint, species)
# Calculate the cumulative sum of n for each correlation
nr_corr_sum <- ddply(nr_corr_sort, c("timepoint", "species"),
                   transform, label_ypos=cumsum(n))


#facet wrap - species
ggplot(data=nr_corr_sum, aes(x=timepoint, y=n, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.6, color="white", size=3.5)+
  scale_y_continuous(limits=c(0, 3000), breaks=c(0, 1000, 2000, 3000))+
  facet_wrap( ~ species) +
  theme(axis.title = element_blank()) +
  labs(title="Number of correlated genes with ZEB2 expression (pearson)", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2")

#facet grid - correlation & species
ggplot(data=nr_corr_sum, aes(x=timepoint, y=n, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= n, label=n), vjust=1.6, color="white", size=3.5)+
  scale_y_continuous(limits=c(0, 2000), breaks=c(0, 1000, 2000))+
  facet_grid(correlation ~ species) +
  theme(axis.title = element_blank()) +
  labs(title="Number of correlated genes with ZEB2 expression (pearson)", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  scale_fill_manual(values=c("#00A86B", "#B90E0A"))

#barplot showing the number of correlation genes throughout the entire time span
nr_corr_all <- data.frame(matrix(, nrow=4, ncol=3))
colnames(nr_corr_all) <- c("correlation", "species","n")
nr_corr_all[,1] <- c(rep("pos", 2), rep("neg",2))
nr_corr_all[,2] <- c(rep(c("human", "gorilla"), 2))
nr_corr_all[,3] <- c(length(rownames(corr_h_ts)[which(sign(corr_h_ts[,1])==1)]),
                     length(rownames(corr_g_ts)[which(sign(corr_g_ts[,1])==1)]),
                     length(rownames(corr_h_ts)[which(sign(corr_h_ts[,1])==-1)]),
                     length(rownames(corr_g_ts)[which(sign(corr_g_ts[,1])==-1)])
                     )
#set the order 
nr_corr_all$species <- factor(nr_corr_all$species, levels = c("human", "gorilla"))
nr_corr_all$correlation <- factor(nr_corr_all$correlation, levels = c("pos", "neg"))

# Sort by timpoint and species
nr_corr_all_sort <- arrange(nr_corr_all, species)
# Calculate the cumulative sum of n for each correlation
nr_corr_all_sum <- ddply(nr_corr_all_sort, c("species"),
                     transform, label_ypos=cumsum(n))
nr_corr_all_sum[,4] <- c(6473, 10232, 3569, 6293)

ggplot(data=nr_corr_all_sum, aes(x=species, y=n, fill=correlation)) +
  geom_bar(stat="identity")  +
  geom_text(aes(y=label_ypos, label=n), vjust=1.6, color="white", size=3.5)+
  theme(axis.title = element_blank()) +
  labs(title="Number of correlated genes with ZEB2 expression \nthroughout the entire time span (pearson)", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  scale_fill_manual(values=c("#00A86B", "#B90E0A"))


#Size of ggplot above: 950x530

#a pie chart showing the number of correlated genes at each time point
nr_corr_total <- aggregate(nr_corr$n, by=list(nr_corr$species,nr_corr$timepoint), FUN=sum)
colnames(nr_corr_total) <- c("species", "timepoint", "n")

#set the order 
nr_corr_total$timepoint <- factor(nr_corr_total$timepoint, levels = c("day0", "day2", "day3", "day5", "day10", "day15", "day25"))

PieDonut(nr_corr_total, aes(species, timepoint, count=n), 
         title = "Number of correlated genes at each time point", titlesize =6,
         showRatioThreshold = 0.001, donutLabelSize = 2.5, r0=0.3,r1=1,r2=1.2,
         showPieName = F)

#a matrix showing the number human-specifically correlated genes at each timepoint without considering pos. and neg correlation
nr_corr_hspc_tp <- data.frame(matrix(, nrow=length(h)*2, ncol=3))
colnames(nr_corr_hspc_tp) <- c("timepoint", "type","n")
nr_corr_hspc_tp[,1] <- paste("day", x, sep = "")
nr_corr_hspc_tp[,2] <- c(rep("h-specific", 7), rep("diff",7))

for (i in 1:nrow(nr_corr_hspc_tp)){
  if (i < 8){
    nr_corr_hspc_tp[i,3] <- sum(!rownames(corr_h[which(!is.na(corr_h[,i])),]) %in% rownames(corr_g[which(!is.na(corr_g[,i])),]))
  } else {
    nr_corr_hspc_tp[i,3] <- sum(!is.na(corr_h[,i-7])) - sum(!rownames(corr_h[which(!is.na(corr_h[,i-7])),]) %in% rownames(corr_g[which(!is.na(corr_g[,i-7])),]))
  }
}

# Sort by timpoint and species
nr_corr_hspc_tp_sort <- arrange(nr_corr_hspc_tp, timepoint)
# Calculate the cumulative sum of n for each correlation
nr_corr_hspc_tp_sum <- ddply(nr_corr_hspc_tp_sort, c("timepoint"),
                     transform, label_ypos=cumsum(n))

#set the order 
nr_corr_hspc_tp_sum$timepoint <- factor(nr_corr_hspc_tp_sum$timepoint,levels = c("day0", "day2", "day3", "day5", "day10", "day15", "day25"))

ggplot(data=nr_corr_hspc_tp_sum, aes(x=timepoint, y=n, fill=type)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.6, color="white", size=3.5)+
  scale_y_continuous(limits=c(0, 3000), breaks=c(0, 1000, 2000,3000))+
  theme(axis.title = element_blank()) +
  labs(title="Number of human-specifically correlated genes with ZEB2 expression (pearson)", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_manual(values=c("#73C2FB", "#8dd35fff"))


#a matrix showing the number human-specifically correlated genes at each timepoint with considering pos. and neg correlation
nr_corr_hspc_tp_pn <- data.frame(matrix(, nrow=length(h)*4, ncol=4))
colnames(nr_corr_hspc_tp_pn) <- c("timepoint", "correlation", "type","n")
nr_corr_hspc_tp_pn[,1] <- paste("day", x, sep = "")
nr_corr_hspc_tp_pn[,2] <- c(rep("pos", 14), rep("neg",14))
nr_corr_hspc_tp_pn[,3] <- c(rep("h-specific", 7), rep("diff",7))

for (i in 1:nrow(nr_corr_hspc_tp_pn)){
  if (i < 15) {
    if (i < 8){
      nr_corr_hspc_tp_pn[i,4] <- sum(!rownames(corr_h)[which(sign(corr_h[,i])==1)] %in% rownames(corr_g)[which(sign(corr_g[,i])==1)])
    } else {
      nr_corr_hspc_tp_pn[i,4] <- as.numeric(length(rownames(corr_h)[which(sign(corr_h[,i-7])==1)])) - sum(!rownames(corr_h)[which(sign(corr_h[,i-7])==1)] %in% rownames(corr_g)[which(sign(corr_g[,i-7])==1)])
    }
  } else {
    if (i < 22) {
      nr_corr_hspc_tp_pn[i,4] <- sum(!rownames(corr_h)[which(sign(corr_h[,i-14])==-1)] %in% rownames(corr_g)[which(sign(corr_g[,i-14])==-1)])
    } else {
      nr_corr_hspc_tp_pn[i,4] <- as.numeric(length(rownames(corr_h)[which(sign(corr_h[,i-21])==-1)])) - sum(!rownames(corr_h)[which(sign(corr_h[,i-21])==-1)] %in% rownames(corr_g)[which(sign(corr_g[,i-21])==-1)])
    }
  }
}

# Sort by timpoint and species
nr_corr_hspc_tp_pn_sort <- arrange(nr_corr_hspc_tp_pn, timepoint, correlation)
# Calculate the cumulative sum of n for each correlation
nr_corr_hspc_tp_pn_sum <- ddply(nr_corr_hspc_tp_pn_sort, c("timepoint", "correlation"),
                             transform, label_ypos=cumsum(n))

#set the order 
nr_corr_hspc_tp_pn_sum$timepoint <- factor(nr_corr_hspc_tp_pn_sum$timepoint,levels = c("day0", "day2", "day3", "day5", "day10", "day15", "day25"))
nr_corr_hspc_tp_pn_sum$correlation <- factor(nr_corr_hspc_tp_pn_sum$correlation,levels = c("pos", "neg"))

ggplot(data=nr_corr_hspc_tp_pn_sum, aes(x=timepoint, y=n, fill=type)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.6, color="white", size=3.5)+
  scale_y_continuous(limits=c(0, 2000), breaks=c(0, 1000, 2000))+
  facet_grid( ~ correlation) +
  theme(axis.title = element_blank()) +
  labs(title="Number of human-specifically correlated genes with ZEB2 expression (pearson)", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_manual(values=c("#73C2FB", "#8dd35fff"))



#a matrix showing the number overlapping genes considering pos & neg at each time point
nr_corr_pn <- data.frame(matrix(, nrow=length(h)*4, ncol=3))
colnames(nr_corr_pn) <- c("timepoint", "correlation","n")
nr_corr_pn[,1] <- paste("day", x, sep = "")
nr_corr_pn[,2] <- c(rep("(+/+)", 7), rep("(+/-)",7), rep("(-/+)", 7), rep("(-/-)", 7))

for (i in 1:nrow(nr_corr_pn)){
  if (i < 15) {
    if (i < 8){
      nr_corr_pn[i,3] <- sum(rownames(corr_h)[which(sign(corr_h[,i+1])==1)] %in% rownames(corr_g)[which(sign(corr_g[,i+1])==1)])
    } else {
      nr_corr_pn[i,3] <- sum(rownames(corr_h)[which(sign(corr_h[,i-6])==1)] %in% rownames(corr_g)[which(sign(corr_g[,i-6])==-1)])
    }
  } else {
    if (i < 22) {
      nr_corr_pn[i,3] <- sum(rownames(corr_h)[which(sign(corr_h[,i-13])==-1)] %in% rownames(corr_g)[which(sign(corr_g[,i-13])==1)])
    } else {
      nr_corr_pn[i,3] <- sum(rownames(corr_h)[which(sign(corr_h[,i-20])==-1)] %in% rownames(corr_g)[which(sign(corr_g[,i-20])==-1)])
    }
  }
}
#set the order 
nr_corr_pn$timepoint <- factor(nr_corr_pn$timepoint,levels = c("day0", "day2", "day3", "day5", "day10", "day15", "day25"))

library(webr)
PieDonut(nr_corr_pn, aes(timepoint, correlation, count=n), 
         title = "Number of overlapping correlated genes with ZEB2", titlesize =6,
         showRatioThreshold = F, donutLabelSize = 2.5, r0=0.3,r1=1,r2=1.2,
         showPieName = F)

#export the list of correlated genes with ZEB2 at each time point.
fullpath <- here("output/correlation_GTF/")
#Human
for (i in 1:ncol(corr_h)) {
  x <- rownames(corr_h)[which(corr_h[,i] < 0)]
  y <- ID_symbol %>% filter(Geneid %in% x)
  write.table(y, file = paste(fullpath,'correlatedGeneswithZEB2_FINAL_', colnames(corr_h)[i],'_neg_n',length(which(corr_h[,i] < 0)),'.tsv', sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
  a <- rownames(corr_h)[which(corr_h[,i] > 0)]
  b <- ID_symbol %>% filter(Geneid %in% a)
  write.table(b, file = paste(fullpath,'correlatedGeneswithZEB2_FINAL_', colnames(corr_h)[i],'_pos_n',length(which(corr_h[,i] > 0)),'.tsv', sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
}

#Gorilla
for (i in 1:ncol(corr_g)) {
  x <- rownames(corr_g)[which(corr_g[,i] < 0)]
  y <- ID_symbol %>% filter(Geneid %in% x)
  write.table(y, file = paste(fullpath,'correlatedGeneswithZEB2_FINAL_', colnames(corr_g)[i],'_neg_n',length(which(corr_g[,i] < 0)),'.tsv', sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
  a <- rownames(corr_g)[which(corr_g[,i] > 0)]
  b <- ID_symbol %>% filter(Geneid %in% a)
  write.table(b, file = paste(fullpath,'correlatedGeneswithZEB2_FINAL_', colnames(corr_g)[i],'_pos_n',length(which(corr_g[,i] > 0)),'.tsv', sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
}



#export the list of correlated genes with ZEB2 thorughout all time
#Human
write.table(ID_symbol %>% filter(Geneid %in% rownames(corr_h_ts)[which(corr_h_ts[,1] < 0)]), 
            file = paste(fullpath,
                            'correlatedGeneswithZEB2_FINAL_H9_', 
                            'alltime_',
                            'neg_n',
                            length(which(corr_h_ts[,1] < 0)),
                            '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")
write.table(ID_symbol %>% filter(Geneid %in% rownames(corr_h_ts)[which(corr_h_ts[,1] > 0)]), 
            file = paste(fullpath,
                         'correlatedGeneswithZEB2_FINAL_H9_', 
                         'alltime_',
                         'pos_n',
                         length(which(corr_h_ts[,1] > 0)),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

#Gorilla
write.table(ID_symbol %>% filter(Geneid %in% rownames(corr_g_ts)[which(corr_g_ts[,1] < 0)]), 
            file = paste(fullpath,
                         'correlatedGeneswithZEB2_FINAL_G1_', 
                         'alltime_',
                         'neg_n',
                         length(which(corr_g_ts[,1] < 0)),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")
write.table(ID_symbol %>% filter(Geneid %in% rownames(corr_g_ts)[which(corr_g_ts[,1] > 0)]), 
            file = paste(fullpath,
                         'correlatedGeneswithZEB2_FINAL_G1_', 
                         'alltime_',
                         'pos_n',
                         length(which(corr_g_ts[,1] > 0)),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

#############################################################################
# ggvenn package for comparison
#############################################################################
library(ggvenn)
library(ggpubr)

#comparison between human & gorilla at each time
List_H <- lapply(1:ncol(corr_h), 
                 function(x){x <- rownames(corr_h)[which(!is.na(corr_h[,x]))]})
names(List_H) <- colnames(corr_h)
List_G <- lapply(1:ncol(corr_g), 
                 function(x){x <- rownames(corr_g)[which(!is.na(corr_g[,x]))]})
names(List_G) <- colnames(corr_g)

#comparison between species at each tp
plotList1 <- lapply(
  1:ncol(corr_h),
  function(CT) {
    x <- ggvenn(
      c(List_H[1], List_G[CT]),
      fill_color = c("#8dd35fff", "#73C2FB"),
      show_percentage = FALSE, set_name_size = 1, text_size = 3
    ) + coord_cartesian(clip="off")
    x
  }
)

#comparison between tp within a species
my_cols <- c('#F8766D','#DB8E00','#94aa01','#03c19f','#00bae3','#629dff','#ff61c3')

plotList1 <- lapply(
  1:ncol(corr_h),
  function(CT) {
    x <- ggvenn(
      c(List_H[1], List_H[CT]),
      fill_color = c(my_cols[1], my_cols[CT]),
      show_percentage = FALSE, set_name_size = 1, text_size = 3
    ) + coord_cartesian(clip="off")
    x
  }
)
 

allplots <- ggarrange(plotlist=c(plotList1, plotList2, 
                                 plotList3, plotList4, 
                                 plotList5, plotList6, plotList7), 
                      ncol = 7, nrow = 7)

allplots

#comparison between human and gorilla with all time correlation
#without pos. and neg. correlation
List_alltime <- list(rownames(corr_h_ts), rownames(corr_g_ts))
names(List_alltime) <- c("human", "gorilla")

ggvenn(
  List_alltime,
  fill_color = c("#8dd35fff", "#73C2FB"),
  show_percentage = TRUE, set_name_size = 3, text_size = 3
)

#with pos. and neg. correlation
List_alltime_pn <- list(rownames(corr_h_ts)[which(sign(corr_h_ts[,1])==1)], 
                        rownames(corr_h_ts)[which(sign(corr_h_ts[,1])==-1)],
                        rownames(corr_g_ts)[which(sign(corr_g_ts[,1])==-1)],
                        rownames(corr_g_ts)[which(sign(corr_g_ts[,1])==1)]
                        )
names(List_alltime_pn) <- c("pos.human", "neg.human", "neg.gorilla", "pos.gorilla")
vennrst <- ggvenn(
  List_alltime_pn,
  fill_color = c("#8dd35fff", "#e0edd8", "#b7dcf7","#73C2FB"),
  show_percentage = TRUE, set_name_size = 4
) + labs(title="Number of correlated genes (+/-)") +
  theme(plot.title = element_text(hjust = 0.5))

List_alltime_pn_cutoff <- list(rownames(corr_h_ts_0.5)[which(sign(corr_h_ts_0.5[,1])==1)], 
                        rownames(corr_h_ts_0.5)[which(sign(corr_h_ts_0.5[,1])==-1)],
                        rownames(corr_g_ts_0.5)[which(sign(corr_g_ts_0.5[,1])==-1)],
                        rownames(corr_g_ts_0.5)[which(sign(corr_g_ts_0.5[,1])==1)]
)

names(List_alltime_pn_cutoff) <- c("pos.human", "neg.human", "neg.gorilla", "pos.gorilla")
vennrst <- ggvenn(
  List_alltime_pn_cutoff,
  fill_color = c("#8dd35fff", "#e0edd8", "#b7dcf7","#73C2FB"),
  show_percentage = TRUE, set_name_size = 4
) + labs(title="Number of correlated genes (+/-)") +
  theme(plot.title = element_text(hjust = 0.5))

path_out <- "output/correlation_GTF/"
#PHNG n = 29
length(intersect(rownames(corr_h_ts_0.5)[which(sign(corr_h_ts_0.5[,1])==1)], rownames(corr_g_ts_0.5)[which(sign(corr_g_ts_0.5[,1])==-1)]))
write.table(ID_symbol %>% filter(Geneid %in% intersect(rownames(corr_h_ts_0.5)[which(sign(corr_h_ts_0.5[,1])==1)], rownames(corr_g_ts_0.5)[which(sign(corr_g_ts_0.5[,1])==-1)])), 
            file = paste(path_out,
                         '2024_08_05_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'PHNG_cutoff0.5_',
                         length(intersect(rownames(corr_h_ts_0.5)[which(sign(corr_h_ts_0.5[,1])==1)], rownames(corr_g_ts_0.5)[which(sign(corr_g_ts_0.5[,1])==-1)])),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

#PHNH n = 31
length(intersect(rownames(corr_h_ts_0.5)[which(sign(corr_h_ts_0.5[,1])==-1)], rownames(corr_g_ts_0.5)[which(sign(corr_g_ts_0.5[,1])==1)]))
write.table(ID_symbol %>% filter(Geneid %in% intersect(rownames(corr_h_ts_0.5)[which(sign(corr_h_ts_0.5[,1])==-1)], rownames(corr_g_ts_0.5)[which(sign(corr_g_ts_0.5[,1])==1)])), 
            file = paste(path_out,
                         '2024_08_05_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'PGNH_cutoff0.5_',
                         length(intersect(rownames(corr_h_ts_0.5)[which(sign(corr_h_ts_0.5[,1])==-1)], rownames(corr_g_ts_0.5)[which(sign(corr_g_ts_0.5[,1])==1)])),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

#############################################################################
#export the gene list for GO analysis 
#############################################################################
###No cutoff
#Human-specific in total n = 7471
write.table(ID_symbol %>% filter(Geneid %in% setdiff(rownames(corr_h_ts), rownames(corr_g_ts))), 
            file = paste(fullpath,
                         '2024_06_28_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'human-specific_',
                         length(setdiff(rownames(corr_h_ts), rownames(corr_g_ts))),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")
#Human-specific pos n = 2802
write.table(ID_symbol %>% filter(Geneid %in% setdiff(rownames(corr_h_ts)[which(corr_h_ts[,1] > 0)], rownames(corr_g_ts))), 
            file = paste(fullpath,
                         '2024_06_28_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'human-specific_onlyPOS_',
                         length(setdiff(rownames(corr_h_ts)[which(corr_h_ts[,1] > 0)], rownames(corr_g_ts))),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")
#Human-specific neg n = 4669
write.table(ID_symbol %>% filter(Geneid %in% setdiff(rownames(corr_h_ts)[which(corr_h_ts[,1] < 0)], rownames(corr_g_ts))), 
            file = paste(fullpath,
                         '2024_06_28_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'human-specific_onlyNEG_',
                         length(setdiff(rownames(corr_h_ts)[which(corr_h_ts[,1] < 0)], rownames(corr_g_ts))),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")
#Overlapping correlated genes with only same direction n = 1727+885 = 2612
write.table(ID_symbol %>% filter(Geneid %in% setdiff(setdiff(intersect(rownames(corr_h_ts), rownames(corr_g_ts)), PH_NG$.), NH_PG$.)), 
            file = paste(fullpath,
                         '2024_06_28_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'overlaps_onlySAMEdirection_',
                         length(setdiff(setdiff(intersect(rownames(corr_h_ts), rownames(corr_g_ts)), PH_NG$.), NH_PG$.)),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

###Cutoff < 0.5
path_out <- "../../Volumes/Extreme SSD/PhD/Drylab/Lancaster/output/correlation_GTF/"
#Human-specific in total n = 5840
length(setdiff(rownames(corr_h_ts_0.5), rownames(corr_g_ts_0.5)))
write.table(ID_symbol %>% filter(Geneid %in% setdiff(rownames(corr_h_ts_0.5), rownames(corr_g_ts_0.5))), 
            file = paste(path_out,
                         '2024_07_23_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'human-specific_cutoff0.5_',
                         length(setdiff(rownames(corr_h_ts_0.5), rownames(corr_g_ts_0.5))),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")
#Human-specific pos n = 2050
length(setdiff(rownames(corr_h_ts_0.5)[which(corr_h_ts_0.5[,1] > 0)], rownames(corr_g_ts_0.5)))
write.table(ID_symbol %>% filter(Geneid %in% setdiff(rownames(corr_h_ts_0.5)[which(corr_h_ts_0.5[,1] > 0)], rownames(corr_g_ts_0.5))), 
            file = paste(path_out,
                         '2024_07_23_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'human-specific_onlyPOS_cutoff0.5_',
                         length(setdiff(rownames(corr_h_ts_0.5)[which(corr_h_ts_0.5[,1] > 0)], rownames(corr_g_ts_0.5))),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")
#Human-specific neg n = 3790
length(setdiff(rownames(corr_h_ts_0.5)[which(corr_h_ts_0.5[,1] < 0)], rownames(corr_g_ts_0.5)))
write.table(ID_symbol %>% filter(Geneid %in% setdiff(rownames(corr_h_ts_0.5)[which(corr_h_ts_0.5[,1] < 0)], rownames(corr_g_ts_0.5))), 
            file = paste(path_out,
                         '2024_07_23_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'human-specific_onlyNEG_cutoff0.5_',
                         length(setdiff(rownames(corr_h_ts_0.5)[which(corr_h_ts_0.5[,1] < 0)], rownames(corr_g_ts_0.5))),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")
#Overlapping correlated genes with only same direction n = 498+1110 = 1608
pp <- intersect(rownames(corr_h_ts_0.5)[which(corr_h_ts_0.5[,1] > 0)],
                rownames(corr_g_ts_0.5)[which(corr_g_ts_0.5[,1] > 0)])
nn <- intersect(rownames(corr_h_ts_0.5)[which(corr_h_ts_0.5[,1] < 0)],
                rownames(corr_g_ts_0.5)[which(corr_g_ts_0.5[,1] < 0)])
length(pp)
length(nn)
write.table(ID_symbol %>% filter(Geneid %in% c(pp, nn)), 
            file = paste(path_out,
                         '2024_07_23_correlatedGeneswithZEB2_', 
                         'alltime_',
                         'overlaps_onlySAMEdirection_cutoff0.5_',
                         length(c(pp, nn)),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")


#############################################################################
# ComplexUpset package for comparison
#############################################################################
library(ComplexUpset)

corr_all_test <- cbind(corr_h, corr_g)
corr_all_test <- as.data.frame(!is.na(corr_all_test)) #41203 rows
corr_all_test <- corr_all_test[which(rowSums(corr_all_test) >= 1),] #17648 rows
upset(corr_all_test, rev(colnames(corr_all_test)), name = "gene", 
      width_ratio=0.1, min_size=50, min_degree=1, sort_sets=FALSE,
      stripes=upset_stripes(
        geom=geom_segment(size=5),
        colors=rev(my_cols)
      ))
upset(corr_all_test, rev(colnames(corr_all_test)), name = "gene", 
      width_ratio=0.1, min_size=3, min_degree=3, sort_sets=FALSE,
      stripes=upset_stripes(
        geom=geom_segment(size=5),
        colors=rev(my_cols)
      ))


#############################################################################
# Gene expression of correlated genes over the entire time window
#############################################################################
#ZEB2 exp: line plot
cts_norm %>% as.data.frame() %>% filter(rownames(cts_norm) == "ENSG00000169554.21") %>%
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>% 
  filter(stringr::str_detect(sample, "H9")) %>%
  mutate(time = tp) %>%
  mutate(time = factor(time, levels = stri_paste("day",x))) %>%
  group_by(time) %>%
  summarise(mean = mean(count)) %>% 
  ggplot(aes(x = time, y = mean, group = 1)) + geom_line() +
  scale_y_log10() + theme_bw() + labs(title="Counts of positively correlated genes (Human)", x="") +
  ylab("normalized count (log)")
  
arrange(corr_h_ts, desc(coef)) %>% slice(2:11) 

#top10 cor in h
top10_h_pos <- arrange(corr_h_ts, desc(coef)) %>% slice(1:11) %>% rownames() #including ZEB2 
top10_h_neg <- arrange(corr_h_ts, coef) %>% slice(c(nrow(corr_h_ts), 1:10)) %>% rownames()
#top10 cor in g
top10_g_pos <- arrange(corr_g_ts, desc(coef)) %>% slice(1:11) %>% rownames() #including ZEB2 
top10_g_neg <- arrange(corr_g_ts, coef) %>% slice(c(nrow(corr_g_ts), 1:10)) %>% rownames()

cts_norm %>% as.data.frame() %>% filter(rownames(cts_norm) == "ENSG00000169554.21") %>%
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count")

cts_norm %>% as.data.frame() %>% slice(match(top10_g_neg, rownames(cts_norm))) %>%
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>% 
  filter(stringr::str_detect(sample, "G1")) %>%
  mutate(time = rep(tp, length(top10_g_neg))) %>%
  mutate(gene = rep(top10_g_neg, each = 21)) %>%
  mutate(time = factor(time, levels = stri_paste("day",x))) %>% 
  group_by(gene, time) %>%
  summarise(mean = mean(count)) %>% 
  mutate(gene = factor(gene, levels = top10_g_neg)) %>%
  ggplot(aes(x = time, y = mean, group = gene)) + geom_line(aes(color=gene)) +
  scale_y_log10() + theme_bw() + labs(title="Counts of negatively correlated genes (Gorilla top10)", x="") +
  ylab("normalized count (log)")

#############################################################################
#with known correlated genes
#############################################################################
#ID2 #ENSG00000115738 #not sig in both
#ZEB2-AS1 #ENSG00000238057.9 #pos sig in both
#CDH1 #ENSG00000039068.19 #neg sig in both
#SNAI1 #ENSG00000124216.4 #neg sig in both
#SNAI2 #ENSG00000019549.13 #pos sig only in human
#NFIL3 #ENSG00000165030.4 #pos sig in both
#BMP4 #ENSG00000125378.16 #neg sig only in human
#ZEB1 #ENSG00000148516.21 #pos only in human
#ZEB1-AS1 #ENSG00000237036.5 #pos only in human
#IRF8 #ENSG00000140968.11 #pos only in gorilla
#GALNT3 #ENSG00000115339.14 #neg sig in both
#ACSL4 #ENSG00000068366 #not sig in both 
#CPT1A #ENSG00000110090 #not sig in both
known_corr <- c("ENSG00000169554.21", #including ZEB2
                "ENSG00000238057.9",
                "ENSG00000039068.19",
                "ENSG00000124216.4",
                "ENSG00000019549.13",
                "ENSG00000165030.4",
                "ENSG00000125378.16",
                "ENSG00000148516.21",
                "ENSG00000237036.5",
                "ENSG00000140968.11",
                "ENSG00000115339.14")

cts_norm %>% as.data.frame() %>% slice(match(known_corr, rownames(cts_norm))) %>%
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>% 
  filter(stringr::str_detect(sample, "G")) %>%
  mutate(time = rep(tp, length(known_corr))) %>%
  mutate(gene = rep(known_corr, each = 21)) %>%
  mutate(time = factor(time, levels = stri_paste("day",x))) %>% 
  group_by(gene, time) %>%
  summarise(mean = mean(count)) %>% 
  mutate(gene = factor(gene, levels = known_corr)) %>%
  ggplot(aes(x = time, y = mean, group = gene)) + geom_line(aes(color=gene)) +
  scale_y_log10() + theme_bw() + labs(title="Counts of known correlated genes (Gorilla)", x="") +
  ylab("normalized count (log)")

#only CDH1 (neg in both), SNAI2 (pos only in h), IRF8 (pos only in g)
OnlyThree <- c("ENSG00000169554.21", #including ZEB2
               "ENSG00000039068.19",
               "ENSG00000019549.13",
               "ENSG00000140968.11")

cts_norm %>% as.data.frame() %>% slice(match(OnlyThree, rownames(cts_norm))) %>%
  pivot_longer(cols = c(1:ncol(cts_norm)), names_to = "sample", values_to = "count") %>% 
  mutate(species = rep(c(rep("Human", 21), rep("Gorilla", 21)), length(OnlyThree))) %>%
  mutate(time = rep(tp, length(OnlyThree)*2)) %>%
  mutate(gene = rep(OnlyThree, each = 21*2)) %>%
  mutate(time = factor(time, levels = stri_paste("day",x))) %>% 
  group_by(gene, time) %>%
  mutate(gene = factor(gene, levels = known_corr)) %>%
  ggplot(aes(x = time, y = count, group = gene)) + geom_point(aes(color=gene)) +
  geom_smooth(aes(color=gene), method=lm , fill="gray", se=TRUE) +
  facet_wrap(species~.) +
  scale_y_log10() + theme_bw() + labs(title="Scatterplot of ZEB2 and its known interactors", x="") +
  ylab("normalized count (log)")
  

#############################################################################
# GO for each time point
#############################################################################
library(clusterProfiler)
library(org.Hs.eg.db)
library(scales)
library(GOSemSim)
library(enrichplot)

go_H0 <- enrichGO(gene=h_day0,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_H2 <- enrichGO(gene=h_day2,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_H3 <- enrichGO(gene=h_day3,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_H5 <- enrichGO(gene=h_day5,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_H10 <- enrichGO(gene=h_day10,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_H15 <- enrichGO(gene=h_day15,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_H25 <- enrichGO(gene=h_day25,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

barplot(go_H0 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_H2 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_H3 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_H5 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_H10 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_H15 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_H25 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')

go_G0 <- enrichGO(gene=g_day0,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_G2 <- enrichGO(gene=g_day2,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_G3 <- enrichGO(gene=g_day3,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_G5 <- enrichGO(gene=g_day5,
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 1)

go_G10 <- enrichGO(gene=g_day10,
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 1)

go_G15 <- enrichGO(gene=g_day15,
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 1)

go_G25 <- enrichGO(gene=g_day25,
                   keyType = "ENSEMBL", 
                   OrgDb = org.Hs.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05, 
                   qvalueCutoff = 1)

barplot(go_G0 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_G2 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_G3 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_G5 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_G10 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_G15 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')
barplot(go_G25 %>% filter(Count > 5), showCategory=30) +scale_y_discrete(labels= label_wrap(50))+ theme(axis.text.y = element_text(size = 8)) + theme(legend.position = 'none')



H0_forGO <- corr_h %>% dplyr::select(1,2) %>% na.omit %>% 
  mutate(gene = rownames(corr_h %>% dplyr::select(1,2) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_h %>% dplyr::select(1,2) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Hu") %>% mutate(day = "0")
colnames(H0_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

H2_forGO <- corr_h %>% dplyr::select(1,3) %>% na.omit %>% 
  mutate(gene = rownames(corr_h %>% dplyr::select(1,3) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_h %>% dplyr::select(1,3) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Hu")%>% mutate(day = "2")
colnames(H2_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

H3_forGO <- corr_h %>% dplyr::select(1,4) %>% na.omit %>% 
  mutate(gene = rownames(corr_h %>% dplyr::select(1,4) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_h %>% dplyr::select(1,4) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Hu") %>% mutate(day = "3")
colnames(H3_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

H5_forGO <- corr_h %>% dplyr::select(1,5) %>% na.omit %>% 
  mutate(gene = rownames(corr_h %>% dplyr::select(1,5) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_h %>% dplyr::select(1,5) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Hu") %>% mutate(day = "5")
colnames(H5_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

H10_forGO <- corr_h %>% dplyr::select(1,6) %>% na.omit %>% 
  mutate(gene = rownames(corr_h %>% dplyr::select(1,6) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_h %>% dplyr::select(1,6) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Hu") %>% mutate(day = "10")
colnames(H10_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

H15_forGO <- corr_h %>% dplyr::select(1,7) %>% na.omit %>% 
  mutate(gene = rownames(corr_h %>% dplyr::select(1,7) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_h %>% dplyr::select(1,7) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Hu") %>% mutate(day = "15")
colnames(H15_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

H25_forGO <- corr_h %>% dplyr::select(1,8) %>% na.omit %>% 
  mutate(gene = rownames(corr_h %>% dplyr::select(1,8) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_h %>% dplyr::select(1,8) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Hu") %>% mutate(day = "25")
colnames(H25_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

G0_forGO <- corr_g %>% dplyr::select(1,2) %>% na.omit %>% 
  mutate(gene = rownames(corr_g %>% dplyr::select(1,2) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_g %>% dplyr::select(1,2) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Go")  %>% mutate(day = "0")
colnames(G0_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

G2_forGO <- corr_g %>% dplyr::select(1,3) %>% na.omit %>% 
  mutate(gene = rownames(corr_g %>% dplyr::select(1,3) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_g %>% dplyr::select(1,3) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Go")  %>% mutate(day = "2")
colnames(G2_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

G3_forGO <- corr_g %>% dplyr::select(1,4) %>% na.omit %>% 
  mutate(gene = rownames(corr_g %>% dplyr::select(1,4) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_g %>% dplyr::select(1,4) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Go")  %>% mutate(day = "3")
colnames(G3_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

G5_forGO <- corr_g %>% dplyr::select(1,5) %>% na.omit %>% 
  mutate(gene = rownames(corr_g %>% dplyr::select(1,5) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_g %>% dplyr::select(1,5) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Go")  %>% mutate(day = "5")
colnames(G5_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

G10_forGO <- corr_g %>% dplyr::select(1,6) %>% na.omit %>% 
  mutate(gene = rownames(corr_g %>% dplyr::select(1,6) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_g %>% dplyr::select(1,6) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Go")  %>% mutate(day = "10")
colnames(G10_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

G15_forGO <- corr_g %>% dplyr::select(1,7) %>% na.omit %>% 
  mutate(gene = rownames(corr_g %>% dplyr::select(1,7) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_g %>% dplyr::select(1,7) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Go")  %>% mutate(day = "15")
colnames(G15_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")

G25_forGO <- corr_g %>% dplyr::select(c(1,8)) %>% na.omit %>% 
  mutate(gene = rownames(corr_g %>% dplyr::select(1,8) %>% na.omit)) %>%
  mutate(reg = ifelse(corr_g %>% dplyr::select(1,8) %>% na.omit %>% dplyr::select(2) > 0, "(+)", "(-)")) %>% mutate(species = "Go")  %>% mutate(day = "25")
colnames(G25_forGO) <- c("symbol" ,"coef", "gene", "reg", "species", "day")


Hu_Go_forGO <- rbind(H0_forGO, H2_forGO, H3_forGO, H5_forGO, H10_forGO, H15_forGO, H25_forGO,
                     G0_forGO, G2_forGO, G3_forGO, G5_forGO, G10_forGO, G15_forGO, G25_forGO)

Hu_forGO <- rbind(H0_forGO, H2_forGO, H3_forGO, H5_forGO, H10_forGO, H15_forGO, H25_forGO)
Go_forGO <- rbind(G0_forGO, G2_forGO, G3_forGO, G5_forGO, G10_forGO, G15_forGO, G25_forGO)


#with SYMBOL
compare.formular_HuGo_sym <- compareCluster(symbol~species+day, 
                                        data =Hu_Go_forGO,
                                        pvalueCutoff = 0.05,
                                        fun = "enrichGO",
                                        OrgDb="org.Hs.eg.db", 
                                        keyType = "SYMBOL", 
                                        ont = "BP")

compare.formular_Hu_sym <- compareCluster(symbol~species+day, 
                                        data =Hu_forGO,
                                        pvalueCutoff = 0.05,
                                        fun = "enrichGO",
                                        OrgDb="org.Hs.eg.db", 
                                        keyType = "SYMBOL", 
                                        ont = "BP")

compare.formular_Go_sym <- compareCluster(symbol~species+day, 
                                        data =Hu_forGO,
                                        pvalueCutoff = 0.05,
                                        fun = "enrichGO",
                                        OrgDb="org.Hs.eg.db", 
                                        keyType = "SYMBOL", 
                                        ont = "BP")

#with ENSEML

compare.formular_HuGo <- compareCluster(gene~species+day, 
                                        data =Hu_Go_forGO,
                                        pvalueCutoff = 0.05,
                                        fun = "enrichGO",
                                        OrgDb="org.Hs.eg.db", 
                                        keyType = "ENSEMBL", 
                                        ont = "BP")

head(compare.formular_HuGo)

compare.formular_Hu <- compareCluster(gene~reg+day, 
                                        data =Hu_forGO,
                                        pvalueCutoff = 0.05,
                                        fun = "enrichGO",
                                        OrgDb="org.Hs.eg.db", 
                                        keyType = "ENSEMBL", 
                                        ont = "BP")
head(compare.formular_Hu)

compare.formular_Hu_day <- compareCluster(gene~day, 
                                      data =Hu_forGO,
                                      pvalueCutoff = 0.05,
                                      fun = "enrichGO",
                                      OrgDb="org.Hs.eg.db", 
                                      keyType = "ENSEMBL", 
                                      ont = "BP")

compare.formular_Go <- compareCluster(gene~reg+day, 
                                      data =Go_forGO,
                                      pvalueCutoff = 0.05,
                                      fun = "enrichGO",
                                      OrgDb="org.Hs.eg.db", 
                                      keyType = "ENSEMBL", 
                                      ont = "BP")

head(compare.formular_Go)

compare.formular_Go_day <- compareCluster(gene~day, 
                                      data =Go_forGO,
                                      pvalueCutoff = 0.05,
                                      fun = "enrichGO",
                                      OrgDb="org.Hs.eg.db", 
                                      keyType = "ENSEMBL", 
                                      ont = "BP")


dot_GO_global <- dotplot(compare.formular_HuGo,font.size=9)
dot_GO_Hu <- dotplot(compare.formular_Hu,font.size=9)
dot_GO_Go <- dotplot(compare.formular_Go,font.size=9)



d <- godata('org.Hs.eg.db', ont="BP")
compare.formular_HuGo <- pairwise_termsim(compare.formular_HuGo, method="Wang", semData = d)
compare.formular_Hu <- pairwise_termsim(compare.formular_Hu, method="Wang", semData = d)
compare.formular_Go <- pairwise_termsim(compare.formular_Go, method="Wang", semData = d)

compare.formular_HuGo_sym <- pairwise_termsim(compare.formular_HuGo_sym, method="Wang", semData = d)
compare.formular_Hu_sym <- pairwise_termsim(compare.formular_Hu_sym, method="Wang", semData = d)
compare.formular_Go_sym <- pairwise_termsim(compare.formular_Go_sym, method="Wang", semData = d)

compare.formular_Hu_day <- pairwise_termsim(compare.formular_Hu_day, method="Wang", semData = d)
compare.formular_Go_day <- pairwise_termsim(compare.formular_Go_sym, method="Wang", semData = d)

emapplot(compare.formular_HuGo)
emapplot(compare.formular_Hu)
emapplot(compare.formular_Go)

emapplot(compare.formular_HuGo_sym)
emapplot(compare.formular_Hu_sym)
emapplot(compare.formular_Go_sym)

emapplot(compare.formular_Hu_day)
emapplot(compare.formular_Go_day)

options(ggrepel.max.overlaps = Inf)
cnetplot(compare.formular_HuGo, node_label = "category")
cnetplot(compare.formular_Hu, node_label = "category")
cnetplot(compare.formular_Go, node_label = "category")

cnetplot(compare.formular_HuGo_sym, node_label = "category")
cnetplot(compare.formular_Hu_sym, node_label = "category")
cnetplot(compare.formular_Go_sym, node_label = "category")

cnetplot(compare.formular_Hu_day, node_label = "category")
cnetplot(compare.formular_Go_day, node_label = "category")


######################################################################
#overlapping GOs
######################################################################
#human (+).day0 - Cadherin binding
#GO:0007156
#GO:0098742
#GO:0016339
######################################################################

######################################################################
#GO ID call function
######################################################################
GOID_call <- function (go_term, sp) 
{
  if(sp == "Hu"){
    go_sp <- c(go_H0, go_H2, go_H3, go_H5,
               go_H10, go_H15, go_H25)
  } else {
    go_sp <- c(go_G0, go_G2, go_G3, go_G5,
               go_G10, go_G15, go_G25)
  }
  
  geneList <- c()
  for (i in 1:length(go_sp)){
    genes <- strsplit(go_sp[[i]]@result[which(go_sp[[i]]@result$ID == go_term),]$geneID, "/")[[1]]
    geneList <- c(geneList, genes)
  }
  print(geneList)
}

######################################################################
#homophilic cell adhesion via plasma membrane adhesion molecules - GO:0007156
######################################################################
Hu_homophilic <- GOID_call(go_term ="GO:0007156", sp ="Hu")
Hu_homophilic <- unique(Hu_homophilic)

Go_homophilic <- GOID_call(go_term ="GO:0007156", sp ="Go")
Go_homophilic <- unique(Go_homophilic)

#overlapping
Hu_homophilic[which(Hu_homophilic %in% Go_homophilic)] #14


plot.homophilic <- lapply(Hu_homophilic[which(Hu_homophilic %in% Go_homophilic)], \(x) {
  y <- plotCounts(dds_MOP, x, c("tp", "sp"), returnData=TRUE)
  y$gene <- x
  return(y)
})

plot.homophilic <- do.call(rbind, plot.homophilic)

ggplot(plot.homophilic, aes(x=tp, y=count,color = sp, group = 1)) +
  geom_point(
    position=position_dodge(width=0.5),
    fill=NA, width=0.25, outlier.shape=NA) +
  stat_summary(aes(group=sp), fun.y=mean, geom="line") +
  facet_wrap(gene~.) +
  scale_y_log10() +
  guides(color=guide_legend("species")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(plot.title = element_text(face="bold", size = 16))

######################################################################
#cell-cell adhesion via plasma-membrane adhesion molecules - GO:0098742
######################################################################
Hu_cc <- GOID_call(go_term ="GO:0098742", sp ="Hu")
Hu_cc <- unique(Hu_cc)

Go_cc <- GOID_call(go_term ="GO:0098742", sp ="Go")
Go_cc <- unique(Go_cc)

#overlapping
Hu_cc[which(Hu_cc %in% Go_cc)] #23


plot.cc <- lapply(Hu_cc[which(Hu_cc %in% Go_cc)], \(x) {
  y <- plotCounts(dds_MOP, x, c("tp", "sp"), returnData=TRUE)
  y$gene <- x
  return(y)
})

plot.cc <- do.call(rbind, plot.cc)

ggplot(plot.cc, aes(x=tp, y=count,color = sp, group = 1)) +
  geom_point(
    position=position_dodge(width=0.5),
    fill=NA, width=0.25, outlier.shape=NA) +
  stat_summary(aes(group=sp), fun.y=mean, geom="line") +
  facet_wrap(gene~.) +
  scale_y_log10() +
  guides(color=guide_legend("species")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(plot.title = element_text(face="bold", size = 16))

######################################################################
#calcium-dependent cell-cell adhesion via plasma membrane cell adhesion molecules - GO:0016339
######################################################################
Hu_kcc <- GOID_call(go_term ="GO:0016339", sp ="Hu")
Hu_kcc <- unique(Hu_kcc)

Go_kcc <- c(
  strsplit(go_G0@result[which(go_G0@result$ID == "GO:0016339"),]$geneID, "/")[[1]],
  strsplit(go_G2@result[which(go_G2@result$ID == "GO:0016339"),]$geneID, "/")[[1]],
  strsplit(go_G3@result[which(go_G3@result$ID == "GO:0016339"),]$geneID, "/")[[1]],
  strsplit(go_G10@result[which(go_G10@result$ID == "GO:0016339"),]$geneID, "/")[[1]],
  strsplit(go_G15@result[which(go_G15@result$ID == "GO:0016339"),]$geneID, "/")[[1]],
  strsplit(go_G25@result[which(go_G25@result$ID == "GO:0016339"),]$geneID, "/")[[1]]
)
Go_kcc <- unique(Go_kcc)

#overlapping
Hu_kcc[which(Hu_kcc %in% Hu_kcc)] 
length(Hu_kcc[which(Hu_kcc %in% Hu_kcc)]) #17

plot.kcc <- lapply(Hu_kcc[which(Hu_kcc %in% Hu_kcc)], \(x) {
  y <- plotCounts(dds_MOP, x, c("tp", "sp"), returnData=TRUE)
  y$gene <- x
  return(y)
})

plot.kcc <- do.call(rbind, plot.kcc)

ggplot(plot.kcc, aes(x=tp, y=count,color = sp, group = 1)) +
  geom_point(
    position=position_dodge(width=0.5),
    fill=NA, width=0.25, outlier.shape=NA) +
  stat_summary(aes(group=sp), fun.y=mean, geom="line") +
  facet_wrap(gene~.) +
  scale_y_log10() +
  guides(color=guide_legend("species")) +
  theme_bw() +
  theme(text = element_text(size = 14)) +
  theme(plot.title = element_text(face="bold", size = 16))
############################################################################################################################################
######################################################################
#Correlation (pre vs post) - general overview
######################################################################
summary_corr <- data.frame(nr = c(length(na.omit(corr_h_peak[,1])), 
                                  length(na.omit(corr_h_peak[,2])),
                                  nrow(na.omit(corr_h_peak)),
                                  length(na.omit(corr_g_peak[,1])),
                                  length(na.omit(corr_g_peak[,2])),
                                  nrow(na.omit(corr_g_peak))),
                              sp = c(rep("human",3), rep("gorilla",3)),
                              pp = c(rep(c("pre","post", "overall"),2))
)

summary_corr %>%
  ggplot(aes(x=factor(pp, levels = c("pre", "post", "overall")), y=as.numeric(nr), fill = factor(sp, levels =c("human", "gorilla")))) +
  geom_bar(stat = "identity",  position=position_dodge()) +
  geom_text(aes(y= nr, label = nr), vjust=1.6, color="white", size=3.5, position = position_dodge(0.9)) +
  labs(title="Number of correlated genes with ZEB2 expression (pearson) pre vs post", x="", y="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  guides(fill=guide_legend(title="species"))

#In case  of correlated genes BEFORE & AFTER ZEB2 peak
corr_peak <- corr_h_peak %>% tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = c(2:3), names_to = "tp", values_to = "coe") %>%
  mutate(species = "human") %>%
  bind_rows(corr_g_peak %>% tibble::rownames_to_column("gene") %>%
              pivot_longer(cols = c(2:3), names_to = "tp", values_to = "coe") %>%
              mutate(species = "gorilla"))

nrow(corr_peak) == (nrow(corr_h_peak)*2)+(nrow(corr_g_peak)*2)

corr_peak %>% group_by(species, tp) %>%
  summarise(non_na_count = sum(!is.na(coe))) %>% 
  ungroup %>%
  add_row(species = "gorilla", tp = "both", non_na_count = corr_g_peak %>% drop_na() %>% nrow(), .before = 3) %>%
  add_row(species = "human", tp = "both", non_na_count = corr_h_peak %>% drop_na() %>% nrow(), .after = 5) %>%
  ggplot(aes(x=factor(tp), y=non_na_count, fill = species)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=non_na_count), vjust=1.6, color="white", position = position_dodge(0.9), size=3.5) %>%
  theme(axis.title = element_blank()) +
  labs(title="Number of correlated genes with ZEB2 expression (pearson) pre vs post", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2")

dir_corr <- function(sp, grp, pos_or_neg){
  if(sp =="Hu"){
    inputdat <- corr_h_peak
  } else {
    inputdat <- corr_g_peak
  }
  unname(table(inputdat[,grp] > 0)[pos_or_neg])
}

summary_corr_pp <- data.frame(nr = c(dir_corr("Hu", 1, 2), dir_corr("Hu", 1, 1),
                                     dir_corr("Hu", 2, 2), dir_corr("Hu", 2, 1),
                                     dir_corr("Go", 1, 2), dir_corr("Go", 1, 1),
                                     dir_corr("Go", 2, 2), dir_corr("Go", 2, 1)),
                              sp = c(rep("human",4), rep("gorilla",4)),
                              pp = c(rep(c(rep("pre",2), rep("post",2)),2)),
                              corr = c(rep(c("pos", "neg"), 4))
)

summary_corr_pp_sum <- ddply(summary_corr_pp, c("pp", "sp"),
                     transform, label_ypos=cumsum(nr))


ggplot(data=summary_corr_pp_sum, aes(x=factor(pp, levels = c("pre", "post")), y=as.numeric(nr), fill=corr)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=nr), vjust=1.6, color="white", size=3.5)+
  scale_y_continuous(limits=c(0, 13000), breaks=c(0, 5000, 10000))+
  facet_wrap( ~ sp) +
  theme(axis.title = element_blank()) +
  labs(title="Number of correlated genes with ZEB2 expression (pearson) \npre vs post", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2")

#export the list of correlated genes 
path_out <- here("output/correlation_GTF/")

#Human
for (i in 1:ncol(corr_h_peak)) {
  x <- rownames(corr_h_peak)[which(corr_h_peak[,i] < 0)]
  y <- ID_symbol %>% filter(Geneid %in% x)
  write.table(y, file = paste("output/correlation_GTF/",'correlatedGeneswithZEB2_peak_H9_', colnames(corr_h_peak)[i],'_neg_n',length(which(corr_h_peak[,i] < 0)),'.txt', sep = ""), row.names = F, col.names = T, quote = F)
  a <- rownames(corr_h_peak)[which(corr_h_peak[,i] > 0)]
  b <- ID_symbol %>% filter(Geneid %in% a)
  write.table(b, file = paste("output/correlation_GTF/",'correlatedGeneswithZEB2_peak_H9_', colnames(corr_h_peak)[i],'_pos_n',length(which(corr_h_peak[,i] > 0)),'.txt', sep = ""), row.names = F, col.names = T, quote = F)
}

#Gorilla
for (i in 1:ncol(corr_g_peak)) {
  x <- rownames(corr_g_peak)[which(corr_g_peak[,i] < 0)]
  y <- ID_symbol %>% filter(Geneid %in% x)
  write.table(y, file = paste("output/correlation_GTF/",'correlatedGeneswithZEB2_peak_G1_', colnames(corr_g_peak)[i],'_neg_n',length(which(corr_g_peak[,i] < 0)),'.txt', sep = ""), row.names = F, col.names = T, quote = F)
  a <- rownames(corr_g_peak)[which(corr_g_peak[,i] > 0)]
  b <- ID_symbol %>% filter(Geneid %in% a)
  write.table(b, file = paste("output/correlation_GTF/",'correlatedGeneswithZEB2_peak_G1_', colnames(corr_g_peak)[i],'_pos_n',length(which(corr_g_peak[,i] > 0)),'.txt', sep = ""), row.names = F, col.names = T, quote = F)
}


######################################################################
#Correlation (pre vs post) - focus on overall correlated genes
######################################################################
#Determining the direction of the correlation (pre vs post)
corr_h_peak_overall <- na.omit(corr_h_peak) %>% 
  mutate(pre = ifelse(pre > 0, "(+)", "(-)")) %>% 
  mutate(post = ifelse(post > 0, "(+)", "(-)"))

corr_g_peak_overall <- na.omit(corr_g_peak) %>% 
  mutate(pre = ifelse(pre > 0, "(+)", "(-)")) %>% 
  mutate(post = ifelse(post > 0, "(+)", "(-)"))

#commends to look into the four variants of the direction changes
#corr_h_peak_overall %>% filter_at(vars(pre, post), all_vars(. %in% "(+)"))  # all pos
#corr_h_peak_overall %>% filter_at(vars(pre, post), all_vars(. %in% "(-)")) # all neg
#corr_h_peak_overall %>% filter_at(vars(pre), all_vars(. %in% "(+)")) %>%
#  filter_at(vars(post), all_vars(. %in% "(-)")) %>% nrow() # pre pos post neg
#corr_h_peak_overall %>% filter_at(vars(pre), all_vars(. %in% "(-)")) %>%
#  filter_at(vars(post), all_vars(. %in% "(+)")) # pre neg post pos

#create a data frame containing species, changes in direction, number of genes
dir_corr <- data.frame(matrix(, nrow = 8, ncol = 3))
colnames(dir_corr) <- c("sp", "ch","n")
dir_corr$sp <- c(rep("Hu",4), rep("Go",4))
dir_corr$ch <- c(rep(c("D/+", "D+", "D/-", "D-"),2))
dir_corr$n <- c(corr_h_peak_overall %>% filter_at(vars(pre, post), all_vars(. %in% "(+)")) %>% nrow(),
                corr_h_peak_overall %>% filter_at(vars(pre), all_vars(. %in% "(-)")) %>%
                  filter_at(vars(post), all_vars(. %in% "(+)")) %>% nrow(),
                corr_h_peak_overall %>% filter_at(vars(pre, post), all_vars(. %in% "(-)")) %>% nrow(),
                corr_h_peak_overall %>% filter_at(vars(pre), all_vars(. %in% "(+)")) %>%
                  filter_at(vars(post), all_vars(. %in% "(-)")) %>% nrow(),
                corr_g_peak_overall %>% filter_at(vars(pre, post), all_vars(. %in% "(+)")) %>% nrow(),
                corr_g_peak_overall %>% filter_at(vars(pre), all_vars(. %in% "(-)")) %>%
                  filter_at(vars(post), all_vars(. %in% "(+)")) %>% nrow(),
                corr_g_peak_overall %>% filter_at(vars(pre, post), all_vars(. %in% "(-)")) %>% nrow(),
                corr_g_peak_overall %>% filter_at(vars(pre), all_vars(. %in% "(+)")) %>%
                  filter_at(vars(post), all_vars(. %in% "(-)")) %>% nrow())

#set the order of the changes
dir_corr$ch <- factor(dir_corr$ch,levels = c("D/+", "D/-", "D+", "D-"))

#pie chart showing the number of overall correlated genes showing the changes in direction of the correlation
PieDonut(dir_corr, aes(sp, ch, count=n), 
         title = "Number of correlated genes including \n changes in the direction of the correlation", titlesize =6,
         showRatioThreshold = F, donutLabelSize = 2.5, r0=0.3,r1=1,r2=1.2,
         showPieName = F)

######################################################################
#Correlation (pre vs post) - GO 1) h vs g, pre vs post, (+) vs (-)
######################################################################
H_pre_pos <- enrichGO(gene=corr_h_peak %>% dplyr::select(pre) %>% na.omit() %>% filter(pre > 0) %>% rownames(),
                  keyType = "ENSEMBL", 
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "none",
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05)

H_pre_neg <- enrichGO(gene=corr_h_peak %>% dplyr::select(pre) %>% na.omit() %>% filter(pre < 0) %>% rownames(),
                       keyType = "ENSEMBL", 
                       OrgDb = org.Hs.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "none",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)

H_post_pos <- enrichGO(gene=corr_h_peak %>% dplyr::select(post) %>% na.omit() %>% filter(post > 0) %>% rownames(),
                       keyType = "ENSEMBL", 
                       OrgDb = org.Hs.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "none",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)

H_post_neg <- enrichGO(gene=corr_h_peak %>% dplyr::select(post) %>% na.omit() %>% filter(post < 0) %>% rownames(),
                       keyType = "ENSEMBL", 
                       OrgDb = org.Hs.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "none",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)

#barplot
barplot(H_pre_pos %>% filter(Count > 5), showCategory=10) +
  scale_y_discrete(labels= label_wrap(50))+ 
  theme(axis.text.y = element_text(size = 8)) + 
  labs(title="GO of correlated genes with ZEB2 expression (H/pre/pos)")
barplot(H_pre_neg %>% filter(Count > 5), showCategory=10) +
  scale_y_discrete(labels= label_wrap(50))+ 
  theme(axis.text.y = element_text(size = 8)) + 
  labs(title="GO of correlated genes with ZEB2 expression (H/pre/neg)")
barplot(H_post_pos %>% filter(Count > 5), showCategory=10) +
  scale_y_discrete(labels= label_wrap(50))+ 
  theme(axis.text.y = element_text(size = 8)) + 
  labs(title="GO of correlated genes with ZEB2 expression (H/post/pos)")
barplot(H_post_neg %>% filter(Count > 5), showCategory=10) +
  scale_y_discrete(labels= label_wrap(50))+ 
  theme(axis.text.y = element_text(size = 8)) + 
  labs(title="GO of correlated genes with ZEB2 expression (H/post/neg)")


G_pre_pos <- enrichGO(gene=corr_g_peak %>% dplyr::select(pre) %>% na.omit() %>% filter(pre > 0) %>% rownames(),
                      keyType = "ENSEMBL", 
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "none",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05)

G_pre_neg <- enrichGO(gene=corr_g_peak %>% dplyr::select(pre) %>% na.omit() %>% filter(pre < 0) %>% rownames(),
                      keyType = "ENSEMBL", 
                      OrgDb = org.Hs.eg.db, 
                      ont = "BP", 
                      pAdjustMethod = "none",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05)

G_post_pos <- enrichGO(gene=corr_g_peak %>% dplyr::select(post) %>% na.omit() %>% filter(post > 0) %>% rownames(),
                       keyType = "ENSEMBL", 
                       OrgDb = org.Hs.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "none",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)

G_post_neg <- enrichGO(gene=corr_g_peak %>% dplyr::select(post) %>% na.omit() %>% filter(post < 0) %>% rownames(),
                       keyType = "ENSEMBL", 
                       OrgDb = org.Hs.eg.db, 
                       ont = "BP", 
                       pAdjustMethod = "none",
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05)

barplot(G_pre_pos %>% filter(Count > 5), showCategory=10) +
  scale_y_discrete(labels= label_wrap(50))+ 
  theme(axis.text.y = element_text(size = 8)) + 
  labs(title="GO of correlated genes with ZEB2 expression (G/pre/pos)")
barplot(G_pre_neg %>% filter(Count > 5), showCategory=10) +
  scale_y_discrete(labels= label_wrap(50))+ 
  theme(axis.text.y = element_text(size = 8)) + 
  labs(title="GO of correlated genes with ZEB2 expression (G/pre/neg)")
barplot(G_post_pos %>% filter(Count > 5), showCategory=10) +
  scale_y_discrete(labels= label_wrap(50))+ 
  theme(axis.text.y = element_text(size = 8)) + 
  labs(title="GO of correlated genes with ZEB2 expression (G/post/pos)")
barplot(G_post_neg %>% filter(Count > 5), showCategory=10) +
  scale_y_discrete(labels= label_wrap(50))+ 
  theme(axis.text.y = element_text(size = 8)) + 
  labs(title="GO of correlated genes with ZEB2 expression (G/post/neg)")

#Heatmap-like functional classification showing gene names 
#import the 123 genes (overlaps between hvsc and hvso) and 79 genes of cluster 3 from b-cell data
g123 <- read.table(here("data/overlap_HvsC_HvsO_123genes.txt"), header = T) #no overlaps AT ALL
g79 <- read.table(here("data/Hspecific_cluster_3_79genes.txt"), header = T) #no overlaps AT ALL
#import the table containing zeb2 target genes
target <- read.csv(here("data/zeb2targetgenes/TFLink_targets_of_O60315.tsv"), sep = "\t" ,header = T)

#Add ensembl id to target table
ENSBL = getBM(attributes = c('external_gene_name', 'ensembl_gene_id'), 
               filters = 'external_gene_name', 
               values = target$Name.Target, 
               mart = ensembl.human)
nrow(ENSBL) #9328 genes
nrow(target) #8662 genes

#sort the gene id in symbol data according to the rownames of data
ENSBL.sort <- ENSBL[match(target$Name.Target, ENSBL$external_gene_name),]
#add a new column of gene symbols sorted by rownames of data
target$Ensembl.Target <- c(ENSBL.sort$ensembl_gene_id)


read_H_pre_pos <- setReadable(H_pre_pos, 'org.Hs.eg.db', 'ENSEMBL')

random.core.genes  <- sapply(
  lapply(
    str_split(read_H_pre_pos@result$geneID, "/"), 
    function(x) intersect(x,target$Name.Target)),
  paste, collapse="/")

#terms annoated to "axonogenesis neuron contraction circulatory"
H_cat <- c("pattern specification process",
           "neuron migration",
           "regulation of neuron projection development",
           "axonogenesis",
           "muscle contraction",
           "vascular process in circulatory system")

# Perform 'hack' by replacing core_enrichment in the gseaResult object.
read_H_pre_pos@result$geneID <- random.core.genes

heatplot(setReadable(H_pre_pos, 'org.Hs.eg.db', 'ENSEMBL'), showCategory="H_cat")
#still too many genes :P


#data prep for merged GO analysis
H_pre_df <- corr_h_peak %>% dplyr::select(coe =pre) %>% na.omit() %>%
  mutate(gene = rownames(corr_h_peak %>% dplyr::select(pre) %>% na.omit())) %>%
  mutate(reg = ifelse(corr_h_peak %>% dplyr::select(pre) %>% na.omit() %>% dplyr::select(1) > 0, "(+)", "(-)")) %>%
  mutate(species = "Hu") %>%
  mutate(tp = "pre")
H_post_df <- corr_h_peak %>% dplyr::select(coe =post) %>% na.omit() %>%
  mutate(gene = rownames(corr_h_peak %>% dplyr::select(post) %>% na.omit())) %>%
  mutate(reg = ifelse(corr_h_peak %>% dplyr::select(post) %>% na.omit() %>% dplyr::select(1) > 0, "(+)", "(-)")) %>%
  mutate(species = "Hu") %>%
  mutate(tp = "post")
G_pre_df <- corr_g_peak %>% dplyr::select(coe =pre) %>% na.omit() %>%
  mutate(gene = rownames(corr_g_peak %>% dplyr::select(pre) %>% na.omit())) %>%
  mutate(reg = ifelse(corr_g_peak %>% dplyr::select(pre) %>% na.omit() %>% dplyr::select(1) > 0, "(+)", "(-)")) %>%
  mutate(species = "Go") %>%
  mutate(tp = "pre")
G_post_df <- corr_g_peak %>% dplyr::select(coe =post) %>% na.omit() %>%
  mutate(gene = rownames(corr_g_peak %>% dplyr::select(post) %>% na.omit())) %>%
  mutate(reg = ifelse(corr_g_peak %>% dplyr::select(post) %>% na.omit() %>% dplyr::select(1) > 0, "(+)", "(-)")) %>%
  mutate(species = "Go") %>%
  mutate(tp = "post")

H_prepost_forGO <- rbind(H_pre_df, H_post_df)
G_prepost_forGO <- rbind(G_pre_df, G_post_df)
HG_prepost_forGO <- rbind(H_pre_df, H_post_df, G_pre_df, G_post_df)

  
#GO run
compare.formular_H_prepost <- compareCluster(gene~tp+reg, 
                                          data =H_prepost_forGO,
                                          pvalueCutoff = 0.05,
                                          fun = "enrichGO",
                                          OrgDb="org.Hs.eg.db", 
                                          keyType = "ENSEMBL", 
                                          ont = "BP")

head(compare.formular_H_prepost)

compare.formular_G_prepost <- compareCluster(gene~tp+reg, 
                                             data =G_prepost_forGO,
                                             pvalueCutoff = 0.05,
                                             fun = "enrichGO",
                                             OrgDb="org.Hs.eg.db", 
                                             keyType = "ENSEMBL", 
                                             ont = "BP")
head(compare.formular_G_prepost)

compare.formular_global_prepost <- compareCluster(gene~tp+reg+species, 
                                             data =HG_prepost_forGO,
                                             pvalueCutoff = 0.05,
                                             fun = "enrichGO",
                                             OrgDb="org.Hs.eg.db", 
                                             keyType = "ENSEMBL", 
                                             ont = "BP")

head(compare.formular_global_prepost)

#dotplot showing GO terms
dot_GO_H_prepost <- dotplot(compare.formular_H_prepost,
                            font.size=9,
                            title = "GO of correlated genes with ZEB2 expression (H)")
dot_GO_G_prepost <- dotplot(compare.formular_G_prepost,
                            font.size=9,
                            title = "GO of correlated genes with ZEB2 expression (G)")
dot_GO_global_prepost <- dotplot(compare.formular_global_prepost,
                                 font.size=9,
                                 title = "GO of correlated genes with ZEB2 expression (H &G)")


d <- godata('org.Hs.eg.db', ont="BP")

compare.formular_Hu_pp <- pairwise_termsim(compare.formular_H_prepost, method="Wang", semData = d)
compare.formular_Go_pp <- pairwise_termsim(compare.formular_G_prepost, method="Wang", semData = d)
compare.formular_HuGo_pp <- pairwise_termsim(compare.formular_global_prepost, method="Wang", semData = d)

emapplot(compare.formular_Hu_pp)
emapplot(compare.formular_Go_pp)
emapplot(compare.formular_HuGo_pp)

options(ggrepel.max.overlaps = Inf)
cnetplot(compare.formular_Hu_pp, node_label = "category")
cnetplot(compare.formular_Go_pp, node_label = "category")
cnetplot(compare.formular_HuGo_pp, node_label = "category")


t1 <- treeplot(compare.formular_Hu_pp, fontsize=2,
               offset.params = list(tiplab = 9, bar_tree = 30))
t2 <- treeplot(compare.formular_Go_pp, fontsize=2,
               offset.params = list(tiplab = 9, bar_tree = 30))
t3 <- treeplot(compare.formular_HuGo_pp, fontsize=2,
               offset.params = list(tiplab = 9, bar_tree = 30))
######################################################################
#Correlation (pre vs post) - GO 2) h vs g - four changes
######################################################################
#corr_h_peak_overall %>% filter_at(vars(pre, post), all_vars(. %in% "(+)"))  # all pos
#corr_h_peak_overall %>% filter_at(vars(pre, post), all_vars(. %in% "(-)")) # all neg
#corr_h_peak_overall %>% filter_at(vars(pre), all_vars(. %in% "(+)")) %>%
#  filter_at(vars(post), all_vars(. %in% "(-)")) %>% nrow() # pre pos post neg
#corr_h_peak_overall %>% filter_at(vars(pre), all_vars(. %in% "(-)")) %>%
#  filter_at(vars(post), all_vars(. %in% "(+)")) # pre neg post pos


#data prep for merged GO analysis
H_Delta <- corr_h_peak_overall %>% 
  mutate(gene = rownames(corr_h_peak_overall)) %>%
  mutate(reg = case_when(corr_h_peak_overall %>% dplyr::select(1) == "(+)" & 
                        corr_h_peak_overall %>% dplyr::select(2) == "(+)" ~ "D/+",
                        corr_h_peak_overall %>% dplyr::select(1) == "(-)" & 
                          corr_h_peak_overall %>% dplyr::select(2) == "(-)" ~ "D/-",
                        corr_h_peak_overall %>% dplyr::select(1) == "(-)" & 
                          corr_h_peak_overall %>% dplyr::select(2) == "(+)" ~ "D+",
                        corr_h_peak_overall %>% dplyr::select(1) == "(+)" & 
                          corr_h_peak_overall %>% dplyr::select(2) == "(-)" ~ "D-")) %>%
  mutate(species = "Hu") 

G_Delta <- corr_g_peak_overall %>% 
  mutate(gene = rownames(corr_g_peak_overall)) %>%
  mutate(reg = case_when(corr_g_peak_overall %>% dplyr::select(1) == "(+)" & 
                           corr_g_peak_overall %>% dplyr::select(2) == "(+)" ~ "D/+",
                         corr_g_peak_overall %>% dplyr::select(1) == "(-)" & 
                           corr_g_peak_overall %>% dplyr::select(2) == "(-)" ~ "D/-",
                         corr_g_peak_overall %>% dplyr::select(1) == "(-)" & 
                           corr_g_peak_overall %>% dplyr::select(2) == "(+)" ~ "D+",
                         corr_g_peak_overall %>% dplyr::select(1) == "(+)" & 
                           corr_g_peak_overall %>% dplyr::select(2) == "(-)" ~ "D-")) %>%
  mutate(species = "Go")

HG_Delta <- rbind(H_Delta, G_Delta)

#GO run
compare.formular_H_Delta <- compareCluster(gene~reg, 
                                             data =H_Delta,
                                             pvalueCutoff = 0.05,
                                             fun = "enrichGO",
                                             OrgDb="org.Hs.eg.db", 
                                             keyType = "ENSEMBL", 
                                             ont = "BP")

compare.formular_G_Delta <- compareCluster(gene~reg, 
                                           data =G_Delta,
                                           pvalueCutoff = 0.05,
                                           fun = "enrichGO",
                                           OrgDb="org.Hs.eg.db", 
                                           keyType = "ENSEMBL", 
                                           ont = "BP")

compare.formular_HG_Delta <- compareCluster(gene~reg+species, 
                                           data =HG_Delta,
                                           pvalueCutoff = 0.05,
                                           fun = "enrichGO",
                                           OrgDb="org.Hs.eg.db", 
                                           keyType = "ENSEMBL", 
                                           ont = "BP")

#dotplot showing GO terms
dot_GO_H_Delta <- dotplot(compare.formular_H_Delta,
                            font.size=9,
                            title = "GO of correlated genes with ZEB2 expression (H, delta)")
dot_GO_G_Delta <- dotplot(compare.formular_G_Delta,
                          font.size=9,
                          title = "GO of correlated genes with ZEB2 expression (G, delta)")
dot_GO_HG_Delta <- dotplot(compare.formular_HG_Delta,
                          font.size=9,
                          title = "GO of correlated genes with ZEB2 expression (HG, delta)")

######################################################################
#positive control with known zeb2 target genes, e.g. CDH1
######################################################################
#e-cadherin gene names and ensembl id
Ecad <- c("CDH1",  "UVO",   "CDHE",  "ECAD",  "LCAM",  "Arc-1", "BCDS1", "CD324", "ENSG00000039068")
#known target genes of zeb2
target_genes <- unique(target$Ensembl.Target)

#correlation test at each time point
corr_h[rownames(corr_h) %in% Ecad,] #negatively correlated only at day 2
corr_g[rownames(corr_g) %in% Ecad,] #NA

corr_h[rownames(corr_h) %in% target_genes,] %>% 
  dplyr::select(-1) %>% filter(if_any(everything(), ~ !is.na(.))) %>% nrow()
#2301 target genes are correlated at least at one time point
corr_g[rownames(corr_g) %in% target_genes,]  %>% 
  dplyr::select(-1) %>% filter(if_any(everything(), ~ !is.na(.))) %>% nrow()
#2526 target genes are correlated at least at one time point

#correlation test before and after zeb2 peak expression
corr_h_peak[rownames(corr_h_peak) %in% Ecad,] #negative correlated in pre
corr_g_peak[rownames(corr_g_peak) %in% Ecad,] #negative correlated in pre

corr_h_peak[rownames(corr_h_peak) %in% target_genes,] %>% filter(if_any(everything(), ~ !is.na(.))) %>% nrow()
#4594 target genes are correlated at least one time range
corr_g_peak[rownames(corr_g_peak) %in% target_genes,]  %>% filter(if_any(everything(), ~ !is.na(.))) %>% nrow()
#4169 target genes are correlated at least one time range

#correlation test before and after zeb2 peak expression - Delta 
H_Delta[rownames(H_Delta) %in% Ecad]
G_Delta[rownames(G_Delta) %in% Ecad]

H_Delta[rownames(H_Delta) %in% target_genes]
H_Delta[rownames(H_Delta) %in% target_genes]

