#############################################################################
#How many lncRNAs in all expressed genes?
#############################################################################
corr_h <- readRDS("output/preliminary/corr_h.rds")
corr_g <- readRDS("output/preliminary/corr_g.rds")

dim(cts_h_filt) #37001    22 - 30061 rows are removed
dim(cts_g_filt) #34826    22 - 32236 rows are removed

sum(gsub("\\..*","",cts_h_filt$Geneid) %in% lncRNA$Ensembl_id) #10015
sum(gsub("\\..*","",cts_g_filt$Geneid) %in% lncRNA$Ensembl_id) #9100


sum(gsub("\\..*","", rownames(corr_h_ts)) %in% lncRNA$Ensembl_id) #2574
sum(gsub("\\..*","", rownames(corr_g_ts)) %in% lncRNA$Ensembl_id) #1437

#############################################################################
#How many lncRNAs in correlated genes at each time point?
#############################################################################
lncRNA #already filtered only for lncRNAs

nr_lncRNAs_tp <- data.frame(matrix(, nrow=length(h)*4, ncol=4))
colnames(nr_lncRNAs_tp) <- c("timepoint", "correlation", "species","n")
nr_lncRNAs_tp[,1] <- paste("day", x, sep = "")
nr_lncRNAs_tp[,2] <- c(rep("pos", 7), rep("neg",7))
nr_lncRNAs_tp[,3] <- c(rep("human", length(h)*2), rep("gorilla", length(h)*2))

for (i in 1:nrow(nr_lncRNAs_tp)){
  if (i < 15) {
    if (i < 8){
      nr_lncRNAs_tp[i,4] <- length(intersect(gsub("\\..*","",rownames(corr_h)[which(corr_h[,i] > 0)]), lncRNA$Ensembl_id))
    } else {
      nr_lncRNAs_tp[i,4] <- length(intersect(gsub("\\..*","",rownames(corr_h)[which(corr_h[,i-7] < 0)]), lncRNA$Ensembl_id))
    }
  } else {
    if (i < 22) {
      nr_lncRNAs_tp[i,4] <- length(intersect(gsub("\\..*","",rownames(corr_g)[which(corr_g[,i-14] > 0)]), lncRNA$Ensembl_id))
    } else {
      nr_lncRNAs_tp[i,4] <- length(intersect(gsub("\\..*","",rownames(corr_g)[which(corr_g[,i-21] < 0)]), lncRNA$Ensembl_id))
    }
  }
}

#change the nr of corr into numeric
nr_lncRNAs_tp$n <- as.numeric(nr_lncRNAs_tp$n)
#set the order 
nr_lncRNAs_tp$timepoint <- factor(nr_lncRNAs_tp$timepoint,levels = c("day0", "day2", "day3", "day5", "day10", "day15", "day25"))
nr_lncRNAs_tp$species <- factor(nr_lncRNAs_tp$species, levels = c("human", "gorilla"))
#add the number of correlated genes 
nr_lncRNAs_tp$totalN <- nr_corr$n

write.table(nr_lncRNAs_tp, file = "output/correlation_GTF/lncRNAsinCorr_eachTP.tsv", 
            col.names = T, row.names = F , quote = F, sep = "\t")

library(plyr)
# Sort by timpoint and species
nr_lncRNAs_tp_sort <- arrange(nr_lncRNAs_tp, timepoint, species)
# Calculate the cumulative sum of n for each correlation
nr_lncRNAs_tp_sum <- ddply(nr_lncRNAs_tp_sort, c("timepoint", "species"),
                     transform, label_ypos=cumsum(n))
# Calculate the percent of lncRNAs
nr_lncRNAs_tp_sum <- nr_lncRNAs_tp_sum %>% mutate(per = n/totalN*100)
# Calculate the cumulative sum of n for each correlation
nr_lncRNAs_tp_sum <- ddply(nr_lncRNAs_tp_sum, c("timepoint", "species"),
                           transform, label_ypos=cumsum(per))
# Average/min/max of percent lncRNAs
mean(nr_lncRNAs_tp_sum$per) #24.2313
min(nr_lncRNAs_tp_sum$per) #14.12587
max(nr_lncRNAs_tp_sum$per) #35.12246

ggplot(data=nr_lncRNAs_tp_sum, aes(x=timepoint, y=n, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=n), vjust=1.6, color="white", size=3.5)+
  scale_y_continuous(limits=c(0, 1100), breaks=c(0, 500, 1000))+
  facet_wrap( ~ species) +
  theme(axis.title = element_blank()) +
  labs(title="Number of lncRNA genes/correlated genes at each time point", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2")

ggplot(data=nr_lncRNAs_tp_sum, aes(x=timepoint, y=per, fill=correlation)) +
  geom_bar(stat="identity") +
  geom_text(aes(y= label_ypos, label=round(per)), vjust=1.6, color="white", size=3.5)+
  scale_y_continuous(limits=c(0, 60), breaks=c(0, 30, 60))+
  facet_wrap( ~ species) +
  theme(axis.title = element_blank()) +
  labs(title="lncRNA genes/correlated genes at each time point", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  theme_bw() +
  labs(y = "%") +
  theme(text = element_text(size = 14))+
  scale_fill_brewer(palette = "Set2")

#############################################################################
#How many lncRNAs in correlated genes all time?
#############################################################################
corr_h_ts

nr_lncRNAs_alltp <- data.frame(matrix(, nrow=4, ncol=3))
colnames(nr_lncRNAs_alltp) <- c("correlation", "species","n")
nr_lncRNAs_alltp[,1] <- c(rep(c("pos", "neg"), 2))
nr_lncRNAs_alltp[,2] <- c(rep("human", 2), rep("gorilla", 2))

for (i in 1:nrow(nr_lncRNAs_alltp)){
  if (i < 3) {
    if (i < 2){
      nr_lncRNAs_alltp[i,3] <- length(intersect(gsub("\\..*","",rownames(corr_h_ts)[which(corr_h_ts[,1] > 0)]), lncRNA$Ensembl_id))
    } else {
      nr_lncRNAs_alltp[i,3] <- length(intersect(gsub("\\..*","",rownames(corr_h_ts)[which(corr_h_ts[,1] < 0)]), lncRNA$Ensembl_id))
    }
  } else {
    if (i < 4) {
      nr_lncRNAs_alltp[i,3] <- length(intersect(gsub("\\..*","",rownames(corr_g_ts)[which(corr_g_ts[,1] > 0)]), lncRNA$Ensembl_id))
    } else {
      nr_lncRNAs_alltp[i,3] <- length(intersect(gsub("\\..*","",rownames(corr_g_ts)[which(corr_g_ts[,1] < 0)]), lncRNA$Ensembl_id))
    }
  }
}

#change the nr of corr into numeric
nr_lncRNAs_alltp$n <- as.numeric(nr_lncRNAs_alltp$n)
#set the order 
nr_lncRNAs_alltp$species <- factor(nr_lncRNAs_alltp$species, levels = c("human", "gorilla"))
#add the number of correlated genes 
nr_lncRNAs_alltp$totalN <- c(length(which(corr_h_ts[,1] > 0)),
                             length(which(corr_h_ts[,1] < 0)),
                             length(which(corr_g_ts[,1] > 0)),
                             length(which(corr_g_ts[,1] < 0)))

write.table(nr_lncRNAs_alltp, file = "output/correlation_GTF/lncRNAsinCorr_allTP.tsv", 
            col.names = T, row.names = F , quote = F, sep = "\t")


#############################################################################
#How many lncRNAs in following three catagories?
#1) Human-specific pos n = 2802
#2) Human-specific neg n = 4669
#3) Overlapping correlated genes with only same direction n = 1727+885 = 2612
#############################################################################

#1) 25 % lncRNAs
sum(gsub("\\..*","", setdiff(rownames(corr_h_ts)[which(corr_h_ts[,1] > 0)], rownames(corr_g_ts))) %in%
  lncRNA$Ensembl_id) #687

#2) 28 % lncRNAs
sum(gsub("\\..*","", setdiff(rownames(corr_h_ts)[which(corr_h_ts[,1] < 0)], rownames(corr_g_ts))) %in%
      lncRNA$Ensembl_id) #1290

#3) 22 % lncRNAs
sum(gsub("\\..*","", setdiff(setdiff(intersect(rownames(corr_h_ts), rownames(corr_g_ts)), PH_NG$.), NH_PG$.)) %in%
      lncRNA$Ensembl_id) #564




