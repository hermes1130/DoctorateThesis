#############################################################################
#correlation test with merged read count table 
#ZEB2 expression vs all other gene in each species

#add gene symbols to ensembl id
#############################################################################
#library calling
library(stats)
library(dplyr)

# Row names are genes, and the first 6 columns are human samples, the last 4 are chimpanzee samples

# Extract human and chimpanzee sample data
human_samples <- cts_norm[, 1:6]
chimp_samples <- cts_norm[, 7:10]

# Define the gene of interest (row name of the gene you want to compare against others)
gene_of_interest <- "ENSG00000169554.23"  # ZEB2

# Function to run Pearson correlation for a set of samples
run_pearson_correlation <- function(data_matrix, target_gene) {
  results <- data.frame(r_value = numeric(nrow(data_matrix)),
                        p_value = numeric(nrow(data_matrix)),
                        row.names = rownames(data_matrix))
  
  target_expression <- data_matrix[target_gene, ]
  
  for (gene in rownames(data_matrix)) {
    if (gene != target_gene) {
      test <- cor.test(target_expression, data_matrix[gene, ], method = "pearson")
      results[gene, "r_value"] <- test$estimate
      results[gene, "p_value"] <- test$p.value
    }
  }
  
  return(results)
}

# Run Pearson correlation for human and chimp samples
human_results <- run_pearson_correlation(human_samples, gene_of_interest)
chimp_results <- run_pearson_correlation(chimp_samples, gene_of_interest)

# Adjust p-values using Benjamini-Hochberg (BH) correction
human_results$adjusted_p_value <- p.adjust(human_results$p_value, method = "BH")
chimp_results$adjusted_p_value <- p.adjust(chimp_results$p_value, method = "BH")

# Filter results based on adjusted p-value cutoff (less than 0.05)
human_filtered <- human_results %>% filter(adjusted_p_value < 0.05) #nothing
chimp_filtered <- chimp_results %>% filter(adjusted_p_value < 0.05) #nothing

# Select only the r_value column as the final result
human_final <- human_filtered %>% select(r_value)
chimp_final <- chimp_filtered %>% select(r_value)

# View the filtered results
print(human_final)
print(chimp_final)


# Create a combined dataframe for both human and chimpanzee samples with a label for the species
human_results$species <- "Human"
chimp_results$species <- "Chimpanzee"

# Combine the data
combined_results <- rbind(human_results, chimp_results)

# Plotting the histogram of adjusted p-values
#png("output/correlation/2024_09_30_Dist_adjustedPval.png", width = 600, height = 800)

ggplot(combined_results, aes(x = adjusted_p_value, fill = species)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Human" = "#8dd35fff", "Chimpanzee" = "#ff7fb0ff")) +
  labs(title = "Distribution of Adjusted p-values",
       x = "Adjusted p-value",
       y = "Frequency") +
  theme_minimal()

#dev.off()

# Plotting the histogram of p-values
#png("output/correlation/2024_09_30_Dist_Pval.png", width = 600, height = 800)

ggplot(combined_results, aes(x = p_value, fill = species)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Human" = "#8dd35fff", "Chimpanzee" = "#ff7fb0ff")) +
  labs(title = "Distribution of p-values",
       x = "p-value",
       y = "Frequency") +
  theme_minimal()

#dev.off()

# Plotting the histogram of r-values
#png("output/correlation/2024_09_30_Dist_Rval.png", width = 600, height = 800)

ggplot(combined_results, aes(x = r_value, fill = species)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Human" = "#8dd35fff", "Chimpanzee" = "#ff7fb0ff")) +
  labs(title = "Distribution of r-values",
       x = "r-value",
       y = "Frequency") +
  theme_minimal()

#dev.off()

# Filter results based on p-value cutoff (less than 0.05)
human_filtered <- human_results %>% filter(p_value < 0.05) 
nrow(human_filtered) #2494
chimp_filtered <- chimp_results %>% filter(p_value < 0.05) 
nrow(chimp_filtered) #3763

# Select only the r_value column as the final result
human_final <- human_filtered %>% select(r_value)
chimp_final <- chimp_filtered %>% select(r_value)

# View the filtered results
print(human_final)
print(chimp_final)

# Paths for the correlation output directory
path_corr <- "output/correlation/"

# Export positively correlated genes for humans
write.table(ID_symbol %>% filter(Geneid %in% rownames(human_final)[which(human_final$r_value > 0)]), 
            file = paste0(path_corr, 'human_pos_corr_pval_n', length(which(human_final$r_value > 0)), '.tsv'), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Export negatively correlated genes for humans
write.table(ID_symbol %>% filter(Geneid %in% rownames(human_final)[which(human_final$r_value < 0)]), 
            file = paste0(path_corr, 'human_neg_corr_pval_n', length(which(human_final$r_value < 0)), '.tsv'), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Export positively correlated genes for chimpanzees
write.table(ID_symbol %>% filter(Geneid %in% rownames(chimp_final)[which(chimp_final$r_value > 0)]), 
            file = paste0(path_corr, 'chimp_pos_corr_pval_n', length(which(chimp_final$r_value > 0)), '.tsv'), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

# Export negatively correlated genes for chimpanzees
write.table(ID_symbol %>% filter(Geneid %in% rownames(chimp_final)[which(chimp_final$r_value < 0)]), 
            file = paste0(path_corr, 'chimp_neg_corr_pval_n', length(which(chimp_final$r_value < 0)), '.tsv'), 
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")



#library calling for venn diagram
library(ggvenn)

#make four lists containing correlation genes
hpos <- rownames(human_final)[which(sign(human_final[,1])==1)]
hneg <- rownames(human_final)[which(sign(human_final[,1])==-1)]
cpos <- rownames(chimp_final)[which(sign(chimp_final[,1])==1)]
cneg <- rownames(chimp_final)[which(sign(chimp_final[,1])==-1)]

#venn diagram 
List <- list(hpos, 
            hneg,
            cneg,
            cneg
)

names(List) <- c("pos.human", "neg.human", "neg.chimp", "pos.chimp")

#png("output/correlation/2024_09_30_ggvenn_nr_corr_pos_neg.png", width = 600, height = 800)
ggvenn(
  List,
  fill_color = c("#8dd35fff", "#e0edd8", "#ead1db","#ff7fb0ff"),
  show_percentage = TRUE, set_name_size = 4
) + labs(title="Number of correlated genes (+/-)") +
  theme(plot.title = element_text(hjust = 0.5))
#dev.off()

#various gene lists from the venn diagram
# Variables that exist only in "x"
hpos_only <- setdiff(hpos, union(cpos, cneg))
length(hpos_only) #1210
hneg_only <- setdiff(hneg, union(cpos, cneg))
length(hneg_only) #957
cpos_only <- setdiff(cpos, union(hpos, hneg))
length(cpos_only) #1617
cneg_only <- setdiff(cneg, union(hpos, hneg))
length(cneg_only) #1819
posInBoth <- intersect(hpos, cpos)
length(posInBoth) #146
negInBoth <- intersect(hneg, cneg)
length(negInBoth) #103
hpos_cneg <- intersect(hpos, cneg)
length(hpos_cneg) #38
hneg_cpos <- intersect(hneg, cpos)
length(hneg_cpos) #39


write.table(ID_symbol %>% filter(Geneid %in% hpos_only), 
            file = paste(path_corr,
                         '2024_10_01_', 
                         'hpos_only_',
                         'pval_n',
                         length(hpos_only),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

write.table(ID_symbol %>% filter(Geneid %in% hneg_only), 
            file = paste(path_corr,
                         '2024_10_01_', 
                         'hneg_only_',
                         'pval_n',
                         length(hneg_only),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

write.table(ID_symbol %>% filter(Geneid %in% cpos_only), 
            file = paste(path_corr,
                         '2024_10_01_', 
                         'cpos_only_',
                         'pval_n',
                         length(cpos_only),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

write.table(ID_symbol %>% filter(Geneid %in% cneg_only), 
            file = paste(path_corr,
                         '2024_10_01_', 
                         'cneg_only_',
                         'pval_n',
                         length(cneg_only),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")   

write.table(ID_symbol %>% filter(Geneid %in% posInBoth), 
            file = paste(path_corr,
                         '2024_10_01_', 
                         'posInBoth_',
                         'pval_n',
                         length(posInBoth),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

write.table(ID_symbol %>% filter(Geneid %in% negInBoth), 
            file = paste(path_corr,
                         '2024_10_01_', 
                         'negInBoth_',
                         'pval_n',
                         length(negInBoth),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

write.table(ID_symbol %>% filter(Geneid %in% hpos_cneg), 
            file = paste(path_corr,
                         '2024_10_01_', 
                         'hpos_cneg_',
                         'pval_n',
                         length(hpos_cneg),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")

write.table(ID_symbol %>% filter(Geneid %in% hneg_cpos), 
            file = paste(path_corr,
                         '2024_10_01_', 
                         'hneg_cpos_',
                         'pval_n',
                         length(hneg_cpos),
                         '.tsv', sep = ""), 
            row.names = F, col.names = T, quote = F, sep = "\t")             

#known ZEB2 target genes
#ID2 #ENSG00000115738.10 ##neg sig only in human
#ZEB2-AS1 #ENSG00000238057.10 #pos sig only in human
#CDH1 #ENSG00000039068.20 
#SNAI1 #ENSG00000124216.4 ##neg sig only in human
#SNAI2 #ENSG00000019549.13 ##neg sig only in human
#NFIL3 #ENSG00000165030.4 
#BMP4 #ENSG00000125378.16 
#ZEB1 #ENSG00000148516.23 
#ZEB1-AS1 #ENSG00000237036.7 
#IRF8 #ENSG00000140968.12 
#GALNT3 #ENSG00000115339.15 
#ACSL4 #ENSG00000068366.21 
#CPT1A #ENSG00000110090.13 
known_corr <- c("ENSG00000169554.23", #including ZEB2
                "ENSG00000115738.10",
                "ENSG00000238057.10",
                "ENSG00000039068.20",
                "ENSG00000124216.4",
                "ENSG00000019549.13",
                "ENSG00000165030.4",
                "ENSG00000125378.16",
                "ENSG00000148516.23",
                "ENSG00000237036.7",
                "ENSG00000140968.12",
                "ENSG00000115339.15",
                "ENSG00000068366.21",
                "ENSG00000110090.13")


#rownames(cts_norm)[grep(paste0("^", "ENSG00000110090"), rownames(cts_norm))]

# Filter genes in human_filtered dataframe where p-value < 0.05
human_significant <- human_filtered[rownames(human_filtered) %in% known_corr & human_filtered$p_value < 0.05, ]

# Filter genes in chimp_filtered dataframe where p-value < 0.05
chimp_significant <- chimp_filtered[rownames(chimp_filtered) %in% known_corr & chimp_filtered$p_value < 0.05, ]

