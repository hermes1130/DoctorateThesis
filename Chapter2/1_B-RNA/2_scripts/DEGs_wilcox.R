#############################################################################
#Alternative to DEseq2
#############################################################################


#############################################################################
#Data prep for Wilcoxon test
#############################################################################
cts_norm <- counts(dds, normalized=TRUE)
dim(cts_norm) #33896    27

cts_norm %>% as.data.frame() %>% filter(rownames(cts_norm) == "ENSG00000169554.21") #ZEB2 is there

#############################################################################
# Data modification - division(ctrl-KD)
#############################################################################
#(+) values = ctrl > KD
#(-) values = ctrl < KD

### create a new data frame to put the divided values
cts_div <- data.frame(matrix(, nrow=nrow(cts_norm), ncol=ncol(cts_norm)))

for(i in 1:ncol(cts_norm)){
  if (i%%3==0){
    cts_div[,i-2] <- cts_norm[,i]-cts_norm[,i-2]
    cts_div[,i-1] <- cts_norm[,i]-cts_norm[,i-1]
  }
}

rownames(cts_div) <- rownames(cts_norm)
colnames(cts_div) <- colnames(cts_norm)

### delete the control column
cts_div <- cts_div %>% dplyr::select(!contains("neg"))
dim(cts_div) #33896    18

#############################################################################
# Pairwise Wilcoxon test
#############################################################################
### prep of result table
rst_test <- data.frame(matrix(, nrow=nrow(cts_div), ncol=4))
colnames(rst_test) <- c("HvsC", "CvsO","empty", "HvsO")
rownames(rst_test) <- rownames(cts_div)

for (i in 1:nrow(cts_div)){
  smaller_than <-
    cts_div[i,] %>% 
    gather(Smp, GeneExp) %>% 
    mutate(Samples = c(rep("Hu", 6), rep("Ch",6), rep("Or", 6))) %>%
    mutate(GeneName = rownames(cts_div[i,])) %>%
    summarize(which(pairwise.wilcox.test(GeneExp,Samples,p.adjust.method = "BH", exact = FALSE)$p.value<0.05))
  result_tbl <- 
    cts_div[i,] %>% 
    gather(Smp, GeneExp) %>% 
    mutate(Samples = c(rep("Hu", 6), rep("Ch",6), rep("Or", 6))) %>%
    mutate(GeneName = rownames(cts_div[i,])) %>%
    summarize(pairwise.wilcox.test(GeneExp,Samples,p.adjust.method = "BH", exact = FALSE)$p.value)
  if (length(smaller_than$`which(...)` >=1)){
    for (j in 1:length(smaller_than$`which(...)`)){
      rst_test[i, smaller_than[j,]] <- result_tbl$`...$p.value`[smaller_than[j,]]
    } 
  } else {
    rm(smaller_than)
    rm(result_tbl)
  }
}

length(na.omit(rst_test$HvsC)) #1406
length(na.omit(rst_test$HvsO)) #1956
length(na.omit(rst_test$CvsO)) #2204

HvsC <- rst_test %>% dplyr::select(HvsC) %>% na.omit
nrow(HvsC) #1406
HvsO <- rst_test %>% dplyr::select(HvsO) %>% na.omit
nrow(HvsO) #1956
CvsO <- rst_test %>% dplyr::select(CvsO) %>% na.omit
nrow(CvsO) #2204


cts_H <- cts_div[, 1:6]
cts_C <- cts_div[, 7:12]
cts_O <- cts_div[, 13:18]
mean_H <- rowMeans(cts_H)
mean_C <- rowMeans(cts_C)
mean_O <- rowMeans(cts_O)

#############################################################################
#Comparision of pairwise wilcoxon test and DE genes detected only in controls
#via DESeq2 (= pure species differences)
#############################################################################
DE_ctrl_HvsC #gene list of default species difference n = 2952
DE_ctrl_HvsO #gene list of default species difference n = 5082

nrow(HvsC) #1406 - pairwise wilcoxon test result
sum(rownames(HvsC) %in% DE_ctrl_HvsC) #223 of 1406 are due to default species difference

nrow(HvsO) #1956 - pairwise wilcoxon test result
sum(rownames(HvsO) %in% DE_ctrl_HvsO) #524 of 1956 are due to default species difference

test <- ID_symbol %>% filter(Geneid %in% DE_ctrl_HvsC)
sum(gsub("\\..*","",test$Geneid) %in% lncRNA$Ensembl_id)


#############################################################################
#table with higher/lower activation comparision
#############################################################################
#add additional information to Hvsx table
#1) activation/repression? a: ctrl-KD >0, r: ctrl-KD <0
#2) higher/lower a/r? higher a/r in human: h>c, lower a/r in human: h<c 
###note! 0 values are assigned to activation
#3) h.c_a.r, h_a.c_r: activated in h repressed in c && h_r.c_a: repressed in h activated in c 
HvsC$a.r_h <- NA
HvsC$a.r_c <- NA
HvsO$a.r_h <- NA
HvsO$a.r_o <- NA

#Exclude genes detected in CvsO
HvsC <- HvsC[which(!rownames(HvsC) %in% rownames(CvsO)),]
nrow(HvsC) #714
HvsO <- HvsO[which(!rownames(HvsO) %in% rownames(CvsO)),]
nrow(HvsO) #901

for (i in rownames(HvsC)){
  ifelse(mean_H[which(names(mean_H)==i)] >= 0, 
         HvsC[which(rownames(HvsC) == i), "a.r_h"] <- "a", 
         HvsC[which(rownames(HvsC) == i), "a.r_h"] <- "r")
}

for (i in rownames(HvsC)){
  ifelse(mean_C[which(names(mean_C)==i)] >= 0, 
         HvsC[which(rownames(HvsC) == i), "a.r_c"] <- "a", 
         HvsC[which(rownames(HvsC) == i), "a.r_c"] <- "r")
}

for (i in rownames(HvsO)){
  ifelse(mean_H[which(names(mean_H)==i)] >= 0, 
         HvsO[which(rownames(HvsO) == i), "a.r_h"] <- "a", 
         HvsO[which(rownames(HvsO) == i), "a.r_h"] <- "r")
}

for (i in rownames(HvsO)){
  ifelse(mean_O[which(names(mean_O)==i)] >= 0, 
         HvsO[which(rownames(HvsO) == i), "a.r_o"] <- "a", 
         HvsO[which(rownames(HvsO) == i), "a.r_o"] <- "r")
}

#New column to add if the gene is activated or repressed by ZEB2
HvsC$h.c_a.r <- NA
HvsO$h.o_a.r <- NA

#Add mean value (KD1 & KD2)
HvsC$mean_H <- mean_H[which(names(mean_H) %in% rownames(HvsC))]
HvsC$mean_C <- mean_C[which(names(mean_C) %in% rownames(HvsC))]

HvsC <- HvsC %>%
  mutate(h.c_a.r = case_when(a.r_h == "a" & a.r_c == "r" ~ "h_a.c_r",
                             a.r_h == "r" & a.r_c == "a" ~ "h_r.c_a"))

#Add mean value (KD1 & KD2)
HvsO$mean_H <- mean_H[which(names(mean_H) %in% rownames(HvsO))]
HvsO$mean_O <- mean_O[which(names(mean_O) %in% rownames(HvsO))]

HvsO <- HvsO %>%
  mutate(h.o_a.r = case_when(a.r_h == "a" & a.r_o == "r" ~ "h_a.o_r",
                             a.r_h == "r" & a.r_o == "a" ~ "h_r.o_a"))

### NA means the gene is either activated or repressed by ZEB2 in both species

#############################################################################
#How many lncRNAs?
#############################################################################
lncRNA <- read.table(here("../data/GRFs_in_GRCh38_complete_list.txt"), header = T)
lncRNA <- lncRNA %>% filter(Gene_type == "lncRNA")

#HvsC
AHRC <- ID_symbol %>% filter(Geneid %in% rownames(HvsC %>% filter(a.r_h == "a" & a.r_c == "r")))
sum(gsub("\\..*","",AHRC$Geneid) %in% lncRNA$Ensembl_id) #69
nrow(AHRC) #282

RHAC <- ID_symbol %>% filter(Geneid %in% rownames(HvsC %>% filter(a.r_h == "r" & a.r_c == "a")))
sum(gsub("\\..*","",RHAC$Geneid) %in% lncRNA$Ensembl_id) #49
nrow(RHAC) #231

AHbiggerAC <- ID_symbol %>% filter(Geneid %in% rownames(HvsC %>% filter(mean_H >= 0, mean_C >= 0) %>% filter(mean_H > mean_C)))
sum(gsub("\\..*","",AHbiggerAC$Geneid) %in% lncRNA$Ensembl_id) #24
nrow(AHbiggerAC) #127

ACbiggerAH <- ID_symbol %>% filter(Geneid %in% rownames(HvsC %>% filter(mean_H >= 0, mean_C >= 0) %>% filter(mean_H < mean_C)))
sum(gsub("\\..*","",ACbiggerAH$Geneid) %in% lncRNA$Ensembl_id) #13
nrow(ACbiggerAH) #42

RCbiggerRH <- ID_symbol %>% filter(Geneid %in% rownames(HvsC %>% filter(mean_H < 0, mean_C < 0) %>% filter(mean_H > mean_C)))
sum(gsub("\\..*","",RCbiggerRH$Geneid) %in% lncRNA$Ensembl_id) #0
nrow(RCbiggerRH) #12

RHbiggerRC <- ID_symbol %>% filter(Geneid %in% rownames(HvsC %>% filter(mean_H < 0, mean_C < 0) %>% filter(mean_H < mean_C)))
sum(gsub("\\..*","",RHbiggerRC$Geneid) %in% lncRNA$Ensembl_id) #4
nrow(RHbiggerRC) #20

#HvsO
AHRO <- ID_symbol %>% filter(Geneid %in% rownames(HvsO %>% filter(a.r_h == "a" & a.r_o == "r")))
sum(gsub("\\..*","",AHRO$Geneid) %in% lncRNA$Ensembl_id) #77
nrow(AHRO) #403

RHAO <- ID_symbol %>% filter(Geneid %in% rownames(HvsO %>% filter(a.r_h == "r" & a.r_o == "a")))
sum(gsub("\\..*","",RHAO$Geneid) %in% lncRNA$Ensembl_id) #45
nrow(RHAO) #297

AHbiggerAO <- ID_symbol %>% filter(Geneid %in% rownames(HvsO %>% filter(mean_H >= 0, mean_O >= 0) %>% filter(mean_H > mean_O)))
sum(gsub("\\..*","",AHbiggerAO$Geneid) %in% lncRNA$Ensembl_id) #36
nrow(AHbiggerAO) #131

AObiggerAH <- ID_symbol %>% filter(Geneid %in% rownames(HvsO %>% filter(mean_H >= 0, mean_O >= 0) %>% filter(mean_H < mean_O)))
sum(gsub("\\..*","",AObiggerAH$Geneid) %in% lncRNA$Ensembl_id) #7
nrow(AObiggerAH) #35

RObiggerRH <- ID_symbol %>% filter(Geneid %in% rownames(HvsO %>% filter(mean_H < 0, mean_O < 0) %>% filter(mean_H > mean_O)))
sum(gsub("\\..*","",RObiggerRH$Geneid) %in% lncRNA$Ensembl_id) #2
nrow(RObiggerRH) #18

RHbiggerRO <- ID_symbol %>% filter(Geneid %in% rownames(HvsO %>% filter(mean_H < 0, mean_O < 0) %>% filter(mean_H < mean_O)))
sum(gsub("\\..*","",RHbiggerRO$Geneid) %in% lncRNA$Ensembl_id) #5
nrow(RHbiggerRO) #17

#These results are based on gene id, but they are different from the results based on gene symbol!!!
#for instance
AHRC <- ID_symbol %>% filter(Geneid %in% rownames(HvsC %>% filter(a.r_h == "a" & a.r_c == "r")))
sum(gsub("\\..*","",AHRC$Geneid) %in% lncRNA$Ensembl_id) #69
nrow(AHRC) #282

AHRC[AHRC$Genesymbol %in% lncRNA$Gene_symbol,] %>% nrow() #14

#############################################################################
#Export gene list in each 6 groups
#############################################################################
#HvsC
#activated in human repressed in chimp
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsC %>% filter(a.r_h == "a" & a.r_c == "r"))),
            file = here("output/results/act_inH_reg_inC_n282.tsv"), 
            quote = F, row.names = F)
#activated in chimp repressed in human
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsC %>% filter(a.r_h == "r" & a.r_c == "a"))),
            file = here("output/results/act_inC_reg_inH_n231.tsv"), 
            quote = F, row.names = F)
#activated in both human & chimp, but human > chimp
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsC %>% filter(mean_H >= 0, mean_C >= 0) %>% filter(mean_H > mean_C))),
            file = here("output/results/highact_inH_vsC_n127.tsv"), 
            quote = F, row.names = F)
#activated in both human & chimp, but chimp > human
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsC %>% filter(mean_H >= 0, mean_C >= 0) %>% filter(mean_H < mean_C))),
            file = here("output/results/highact_inC_vsH_n42.tsv"), 
            quote = F, row.names = F)
#repressed in both human & chimp, but human > chimp
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsC %>% filter(mean_H < 0, mean_C < 0) %>% filter(mean_H < mean_C))),
            file = here("output/results/highreg_inH_vsC_n12.tsv"), 
            quote = F, row.names = F)
#repressed in both human & chimp, but chimp > human
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsC %>% filter(mean_H < 0, mean_C < 0) %>% filter(mean_H > mean_C))),
            file = here("output/results/highreg_inC_vsH_n20.tsv"), 
            quote = F, row.names = F)

#HvsO

#activated in human repressed in orang
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsO %>% filter(a.r_h == "a" & a.r_o == "r"))),
            file = here("output/results/act_inH_reg_inO_n403.tsv"), 
            quote = F, row.names = F)
#activated in orang repressed in human
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsO %>% filter(a.r_h == "r" & a.r_o == "a"))),
            file = here("output/results/act_inO_reg_inH_n297.tsv"), 
            quote = F, row.names = F)
#activated in both human & orang, but human > orang
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsO %>% filter(mean_H >= 0, mean_O >= 0) %>% filter(mean_H > mean_O))),
            file = here("output/results/highact_inH_vsO_n131.tsv"), 
            quote = F, row.names = F)
#activated in both human & orang, but orang > human
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsO %>% filter(mean_H >= 0, mean_O >= 0) %>% filter(mean_H < mean_O))),
            file = here("output/results/highact_inO_vsH_n35.tsv"), 
            quote = F, row.names = F)
#repressed in both human & orang, but human > orang
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsO %>% filter(mean_H < 0, mean_O < 0) %>% filter(mean_H < mean_O))),
            file = here("output/results/highreg_inH_vsO_n17.tsv"), 
            quote = F, row.names = F)
#repressed in both human & orang, but orang > human
write.table(ID_symbol %>% 
              filter(Geneid %in% rownames(HvsO %>% filter(mean_H < 0, mean_O < 0) %>% filter(mean_H > mean_O))),
            file = here("output/results/highreg_inO_vsH_n18.tsv"), 
            quote = F, row.names = F)

#human-specific n = 417
write.table(ID_symbol %>% 
              filter(Geneid %in% setdiff(rownames(HvsC)[which(rownames(HvsC) %in% rownames(HvsO))], rownames(CvsO))),
            file = here("output/results/hspecifc_n417.tsv"), 
            quote = F, row.names = F)










