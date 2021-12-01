# ----
# Look at the pairwise correlation between HIF target gene list genes
# ----

rm(list=ls())
library(tidyverse)
library(ComplexHeatmap)

# HIF target gene list 
# removed endothelial markers (FLT1 and DLL4),
# removed other markers that were missing from the Pseudbulk data (likely removed for having low expression across all cells)
HIF_genes <- c("SLC2A1",
               "PLIN2",
               "CA12",
               "FLG", #
               "IL6",
               "ADM",
               "VEGFA",
               "BNIP3",
               "HK1",
               "HK2",
               "PFKM",
               "ALDOA",
               "PGK1",
               "LDHA",
               "NOS2",
               "ABL2",
               "EPO", #
               "POU5F1", #
               "SCGB3A1", #
               "TGFA",
               "CCND1",
               "DLL4", #
               "ANGPT2",
               "EDN1",
               "EGLN3",
               "EGLN1",
               "TF",
               "TGFB3",
               "FLT1", #
               "ARNT", 
               "EPAS1")

# this dataframe was originally generated in the dotplot_chrom_diff_genes.R script
# it indicates fraction of cells with > 1 read count for each gene

fraction_cells_expressing <- read_tsv("Data/pseudobulk_fraction_cells_expressing.tsv")
# a column also indicates whether or not the gene is expressed in >= 20% of cells 

fraction_cells_expressing_HIF <- fraction_cells_expressing %>%
  dplyr::select(-gene_over_20pct_in_at_least_one_sample) %>% # remove the old version of the column 
  # get genes of interest
  filter(gene %in% HIF_genes) %>% 
  # make column describing which genes are expressed in at least 20% of tumour cells in at least 1 sample
  rowwise() %>%
  mutate(gene_20pct_or_over_in_at_least_one_sample = if_else(sum(c_across(`E035_Tumour`:`E007_Tumour`) >= 0.2) > 0, TRUE, FALSE)) %>% 
  ungroup()

# filter out the HIF genes that are below the expression threshold
HIF_genes_tumour <- fraction_cells_expressing_HIF %>%
  filter(gene_20pct_or_over_in_at_least_one_sample == TRUE) %>%
  pull(gene)

# ----
# HIF target correlation in PCPG
# ----

# read in the pseudobulk data 
pseudobulk_pcpg <- read_csv("Data/Pseudobulk_log2_CPM.csv")
pseudobulk_pcpg <- pseudobulk_pcpg %>%
  dplyr::rename("Gene" = "gene") %>%
  dplyr::select(Gene, ends_with("Chromaffin.cells")) %>%  # get just the chromaffin cells from each tumour and normal sample 
  data.frame() 
# remove the cell type label from the sample name as now they are all chromaffin / tumour 
colnames(pseudobulk_pcpg) <- gsub(pattern = "_Chromaffin.cells",
                                  replacement = "",
                                  colnames(pseudobulk_pcpg)) 

# remove the normal adrenal samples E240 and E243
pseudobulk_pcpg <- pseudobulk_pcpg[,!colnames(pseudobulk_pcpg) %in% c("E240", "E243")]

rownames(pseudobulk_pcpg) <- pseudobulk_pcpg$Gene
pseudobulk_pcpg$Gene <- NULL
pseudobulk_pcpg <- as.matrix(pseudobulk_pcpg)
HIF_gene_expression <- pseudobulk_pcpg[HIF_genes_tumour, ]

# calculate pairwise correlation for each HIF target gene 
# Spearman correlation 
spearman_correlations <-  list()
for (i in 1:nrow(HIF_gene_expression)){
  gene <- HIF_gene_expression[i,]
  spearman_gene <- cor(gene, t(HIF_gene_expression),
                       method = "pearson")
  spearman_correlations[[rownames(HIF_gene_expression)[i]]] <- spearman_gene
}
spearman_correlations_mat <- Reduce(f = rbind, x = spearman_correlations)
rownames(spearman_correlations_mat) <- rownames(HIF_gene_expression)

# write.csv(spearman_correlations_mat,
#           "Results/VEGFA_spearman_correlation_HIF_target_genes.tsv")
