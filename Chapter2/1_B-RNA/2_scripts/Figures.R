library("ggvenn")
library("pheatmap")
library("tibble")
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")


#############################################################################
#venn diagram
#############################################################################
#make a list with main effect
maineff <- list(
  Human =rownames(Main_H_sig),
  Chimp =rownames(Main_C_sig),
  Orang =rownames(Main_O_sig)
)


#make a list with main effect difference
maineffdiff <- list(
  HuvsCh =rownames(condDiff_HvsC_sig),
  HuvsOr =rownames(condDiff_HvsO_sig),
  ChvsOr =rownames(condDiff_CvsO_sig)
)
capture.output(sig_genes0.01, file = "sig_genes0.01.csv")
#plotting
ggvenn(
  maineff, 
  fill_color = c("#8dd35fff","#ff80b2ff", "#ff9955ff"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE) +
  ggtitle("DEGs in KD compared to ctrl") + 
  theme(plot.title = element_text(size=18,  face="bold",hjust = 0.5))

ggvenn(
  maineffdiff, 
  fill_color = c("#8dd35fff", "#ff80b2ff", "#ff9955ff"),
  stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE) +
  ggtitle("DEGs in NHP compared to human") + 
  theme(plot.title = element_text(size=18,  face="bold",hjust = 0.5))

cowplot::plot_grid(v0.05, v0.01, ncol=2, labels=LETTERS[1:2])

#############################################################################
#Clustering + heatmap - with pos. & neg. regulation
#############################################################################
### DEGs adj p.val < 0.05
sig0.01 <- c(rownames(rst_test %>% filter(HvsC < 0.01)), 
             rownames(rst_test %>% filter(HvsO < 0.01)), 
             rownames(rst_test %>% filter(CvsO < 0.01)))
length(sig0.01) #866
### remove duplicated genes
sig0.01 <-unique(sig0.01) 
length(sig0.01) #494
### DEGs adj p.val < 0.05
sig0.05 <- c(rownames(rst_test[!is.na(rst_test$HvsC),]), 
             rownames(rst_test[!is.na(rst_test$HvsO),]), 
             rownames(rst_test[!is.na(rst_test$CvsO),]))
length(sig0.05) #5566
### remove duplicated genes
sig0.05 <-unique(sig0.05) 
length(sig0.05) #3402

### read count table with genes detected via Wilcoxon adj p.val < 0.01 & < 0.05
cts_sig0.01 <- cts_div[rownames(cts_div) %in% sig0.01, ] #%>% filter(!Orang_JingJing_Zeb2_3 > 8000)
cts_sig0.05 <- cts_div[rownames(cts_div) %in% sig0.05, ] #%>% filter(!Orang_JingJing_Zeb2_3 > 8000)

### histogram before normalization 
hist(as.matrix(cts_sig0.01))
hist(as.matrix(cts_sig0.05))

range(cts_sig0.01) #-3908.976 - 65708.19
range(cts_sig0.05) #-62601.24 - 143072.2

### To transform it to a log scale, 
### first find the log of the positive number then multiply it by its sign
transform_to_log_scale <- function(x){
  ifelse(x==0, log(1), (sign(x)) * (log2(abs(x)+1)))
}

### log normalization
cts_sig0.01_norm <- sapply(cts_sig0.01, transform_to_log_scale)
hist(cts_sig0.01_norm)
range(cts_sig0.01_norm) #-11.93294 - 16.00381

cts_sig0.05_norm <- sapply(cts_sig0.05, transform_to_log_scale)
hist(cts_sig0.05_norm)
range(cts_sig0.05_norm) #-15.93393 - 17.12639

### create a matrix
rownames(cts_sig0.01_norm) <- rownames(cts_sig0.01)
rownames(cts_sig0.05_norm) <- rownames(cts_sig0.05)
hclust_matrix0.01 <- as.matrix(cts_sig0.01_norm)
hclust_matrix0.05 <- as.matrix(cts_sig0.05_norm)

### scale the data (to obtain z-scores)
#it will scale each column by its mean and standard deviation. 
#However, we want to scale the expression of our genes, which are the rows of the matrix! 
#So, we need to do a little gymnastics here
#first transpose our matrix, then scale, then transpose it back again.
hclust_matrix0.01 <- hclust_matrix0.01 %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()


hclust_matrix0.05 <- hclust_matrix0.05 %>% 
  # transpose the matrix so genes are as columns
  t() %>% 
  # apply scalling to each column of the matrix (genes)
  scale() %>% 
  # transpose back so genes are as rows again
  t()

### calculate the distance between each gene (row) in our matrix
gene_dist0.01 <- dist(hclust_matrix0.01)
gene_dist0.05 <- dist(hclust_matrix0.05)

### hierarchical clustering using hclust()
gene_hclust0.01 <- hclust(gene_dist0.01, method = "complete")
gene_hclust0.05 <- hclust(gene_dist0.05, method = "complete")

#heatmap
#labels <- unique(trans_cts_mean0.01 %>% filter(cluster == "1") %>% pull(gene))
pheatmap(hclust_matrix0.01, cluster_rows=TRUE,  show_rownames=FALSE, treeheight_row = 100,
         cluster_cols=FALSE, show_colnames = FALSE)+ theme(legend.position = 'none') 

pheatmap(hclust_matrix0.05, cluster_rows=TRUE, show_rownames=FALSE, treeheight_row = 100,
         cluster_cols=FALSE, show_colnames = FALSE)+ theme(legend.position = 'none') 



# The default `plot()` function can be used to produce a simple dendrogram
plot(gene_hclust0.01, labels = FALSE)
abline(h = 8, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

plot(gene_hclust0.05, labels = FALSE)
abline(h = 8, col = "brown", lwd = 2) # add horizontal line to illustrate cutting dendrogram

#dendrogram “cutting”
gene_cluster0.01 <- cutree(gene_hclust0.01, k = 6) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)

gene_cluster0.05 <- cutree(gene_hclust0.05, k = 6) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)

head(gene_cluster0.01)
head(gene_cluster0.05)


sample <- c("Human_R1_zeb2_2", "Human_R1_zeb2_3", 
            "Human_R2_zeb2_2", "Human_R2_zeb2_3", 
            "Human_R3_zeb2_2", "Human_R3_zeb2_3", 
            "Chimp_Judith_Zeb2_2", "Chimp_Judith_Zeb2_3", 
            "Chimp_Leo_Zeb2_2", "Chimp_Leo_Zeb2_3", 
            "Chimp_Maryke_Zeb2_2", "Chimp_Maryke_Zeb2_3", 
            "Orang_Guchi_Zeb2_2", "Orang_Guchi_Zeb2_3", 
            "Orang_Jaqo_Zeb2_2", "Orang_Jaqo_Zeb2_3", 
            "Orang_JingJing_Zeb2_2", "Orang_JingJing_Zeb2_3")
condition <- c(rep(c("KD_1", "KD_2"), 9))
species <- c(rep("Hu", 6), 
             rep("Ch", 6),
             rep("Or", 6))
replicat <- c(rep(c(rep("A",2),
                    rep("B",2),
                    rep("C",2)),3))
sample_info <- data.frame(sample, condition, species, replicat)

### summarise counts 
### add a new column of genes
cts_sig0.01_norm <- cbind(gene = rownames(cts_sig0.01_norm), cts_sig0.01_norm)
cts_sig0.05_norm <- cbind(gene = rownames(cts_sig0.05_norm), cts_sig0.05_norm)

### log normalized count table -> scale -> mean per each replicate
trans_cts_mean0.01 <- as.data.frame(cts_sig0.01_norm) %>% 
  # convert to long format
  pivot_longer(cols = c(2:19), names_to = "sample", values_to = "cts")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample")) %>%
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(cts_scaled = (as.numeric(cts) - mean(as.numeric(cts)))/sd(as.numeric(cts))) %>% 
  # for each gene, species and replicat
  group_by(gene, species, replicat) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled)) %>% 
  ungroup()

trans_cts_mean0.05 <- as.data.frame(cts_sig0.05_norm) %>% 
  # convert to long format
  pivot_longer(cols = c(2:19), names_to = "sample", values_to = "cts")  %>% 
  # join with sample info table
  full_join(sample_info, by = ("sample")) %>%
  # for each gene
  group_by(gene) %>% 
  # scale the cts column
  mutate(cts_scaled = (as.numeric(cts) - mean(as.numeric(cts)))/sd(as.numeric(cts))) %>% 
  # for each gene, species and replicat
  group_by(gene, species, replicat) %>%
  # calculate the mean (scaled) cts
  summarise(mean_cts_scaled = mean(cts_scaled)) %>% 
  ungroup()


#join our table with gene clusters to the table of gene summarised counts from above  
trans_cts_mean0.01 <- trans_cts_mean0.01 %>% 
  inner_join(gene_cluster0.01, by = "gene")

head(trans_cts_mean0.01)
trans_cts_mean0.01$species <- factor(trans_cts_mean0.01$species, levels=c("Hu", "Ch", "Or"))

trans_cts_mean0.05 <- trans_cts_mean0.05 %>% 
  inner_join(gene_cluster0.05, by = "gene")

head(trans_cts_mean0.05)
trans_cts_mean0.05$species <- factor(trans_cts_mean0.05$species, levels=c("Hu", "Ch", "Or"))

#plotting
cl_0.01 <- trans_cts_mean0.01 %>% 
  ggplot(aes(species, mean_cts_scaled)) +
  geom_line(aes(group = gene), alpha = 0.3) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1.5, 
            aes(group = 1)) +
  labs(title="Gene expression in each cluster (DEGs adj p.val < 0.01)", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  ylab("Normlized gene expression") +
  theme_bw() +
  facet_grid(cols = vars(cluster))

cl_0.05 <- trans_cts_mean0.05 %>% 
  ggplot(aes(species, mean_cts_scaled)) +
  geom_line(aes(group = gene), alpha = 0.3) +
  geom_line(stat = "summary", fun = "median", colour = "brown", size = 1.5, 
            aes(group = 1)) +
  labs(title="Gene expression in each cluster (DEGs adj p.val < 0.05)", x="") +
  theme(plot.title = element_text(face="bold", size = 16)) +
  ylab("Normlized gene expression") +
  theme_bw() +
  facet_grid(cols = vars(cluster))


#heatmap again for the final version including clusters
gene_col0.05 <- data.frame(gene_cluster0.05)
rownames(gene_col0.05) <- gene_cluster0.05$gene
gene_col0.05 <- gene_col0.05 %>% dplyr::select(cluster)

pheatmap(hclust_matrix0.05, cluster_rows=TRUE, show_rownames=FALSE, 
         cluster_cols=TRUE, show_colnames = TRUE, cutree_rows = 6,
         annotation_row = gene_col0.05)







