rm(list=ls())

library(tidyverse)
library(Seurat)
library(circlize)
library(RColorBrewer)
library(ComplexHeatmap)
source("Figures/dotplot_functions.R")
source("singlecell_colour_palettes.R")

# Read in bulk batch normalised expression
bulk_norm_expression <- read_tsv("Data/bulk_expression_batch_normalised.tsv")

# Read in the bulk metadata
bulk_md <- read_csv("Data/Table S3. Pheo-atlas metadata.csv") %>%
  mutate(Sample = colnames(bulk_norm_expression[,2:ncol(bulk_norm_expression)]))

# check samples are in the same order
table(colnames(bulk_norm_expression[,2:ncol(bulk_norm_expression)]) == bulk_md$Sample)

# read in tumour vs rest DGE
de_tables <- read_tsv("/data/gpfs/projects/punim0010/projects/Pattison_projects/PPGL_sc_RNA_seq_compendium/Bulk/Bulk DE/bulk_de_by_subtype_reversed.tsv") %>%
  filter(!Cluster %in% c("Cortical.admixture")) %>%
  # use logFC threshold of +- 0.5 and pval of 0.05
  filter(adj.P.Val < 0.05) %>%
  filter(abs(logFC) > 0.5) %>%
  # make column describing if over or underexpressed 
  mutate(direction = ifelse(logFC > 0 , "Up", "Down")) %>% 
  mutate(gene_direction = paste(gene, direction, sep = "_")) %>% 
  group_by(gene) %>%
  mutate(count = n()) %>%
  ungroup() %>% 
  # rename the clusters so that they are consistent with the bulk metadata
  mutate(Cluster = recode(Cluster,
    "SDHx_HN" = "SDHx (H&N)",
    "PH_NOS" = "PH-NOS"))

# read in tumour vs normal DGE 
tumour_normal <- read_csv("/data/gpfs/projects/punim0010/projects/Pattison_projects/PPGL_sc_RNA_seq_compendium/Tables/Paper_tables/Pseudobulk tumour vs chromaffin DGE.csv") %>%
  filter(Contrast %in% c("C1Ai","C1Aii","C1Bi", "C1Bii","C2A", "C2Bi", "C2Bii")) %>%
  filter(adj.P.Val < 0.05) %>%
  filter(abs(logFC) > 0.5) %>%
  # make column describing if over or underexpressed 
  mutate(direction = ifelse(logFC > 0 , "Up", "Down")) %>%
  mutate(gene_direction = paste(Gene, direction, sep = "_"))

# make a vector containing old and new cluster names 
names_df <- data.frame(Cluster = c(all.clusters, "Multiple"),
                       Comete = c(all.clusters2, "Multiple"))

DE_in_pseudobulk_tumour_v_normal.genes <- tumour_normal$gene_direction

# Get the genes that are DE in only one sample
one_sample <- de_tables %>%
  ungroup() %>%
  mutate(DE_in_pseudobulk = ifelse(gene_direction %in% tumour_normal$gene_direction, "Yes", "No")) %>%
  mutate(DE_in_pseudobulk = factor(DE_in_pseudobulk, levels = c("Yes", "No"))) %>%
  left_join(names_df) %>%
  mutate(Comete = replace(Comete, count > 1, "Multiple")) %>%
  arrange(DE_in_pseudobulk) %>%
  filter(!duplicated(gene))

# Get genes unique to a subtype (bulk) and DE in both bulk and PB
DE_in_both <- filter(one_sample, 
                     DE_in_pseudobulk == "Yes")

# Genes that will end up in the heatmap
Heatmap_genes <- de_tables %>%
  filter(gene %in% one_sample$gene) %>%
  # Join on the mean proportion of cells expressing this gene in each cell type 
  mutate(Gene_is_DE_in_tumour_PB = ifelse(gene %in% DE_in_both$gene, "Yes", "No"))

subgroups.lvl <- c("C1Ai", "C1Aii", "C1Bi", "C1Bii", "C2A", "C2Bi", "C2Bii", "Multiple")

bulk_md <- bulk_md %>%
  mutate(Subtype_old = Subtype) %>% 
  rename(Cluster = Subtype)

anno <- bulk_md %>%
  left_join(names_df, by = "Cluster") %>%
  rename(Subtype = Comete) %>%
  filter(!Subtype %in% c("C2C", "Normal")) %>%
  mutate(Subtype = recode(Subtype, !!!subtype_key1)) %>% 
  mutate(Subtype = factor(Subtype, levels = subtype_key1)) %>%
  arrange(Subtype)

# Convert to matrix, select genes and z score scale
bulk_norm_expression_mat <- as.matrix(bulk_norm_expression[,2:ncol(bulk_norm_expression)])
rownames(bulk_norm_expression_mat) <- bulk_norm_expression$Gene
# Select columns and rows
bulk_norm_expression_mat <- bulk_norm_expression_mat[one_sample$gene,anno$Sample]
# Z score transform
bulk_norm_expression_mat <- t(scale(t(bulk_norm_expression_mat)))

# Set NAs to 0
bulk_norm_expression_mat[is.na(bulk_norm_expression_mat)] <- 0

# Top annotation describing the subtype
gen_anno1 <- data.frame(Subtype = anno$Subtype)
top_ano <- HeatmapAnnotation(df = gen_anno1,
                             col = list("Subtype" = c(subtype_genotype_cols[2:8],
                                                      "Multiple" = "lightgrey")),
                             show_annotation_name = FALSE)

# Left annotation
gen_anno <- data.frame(subtype_specific = recode(one_sample$Comete, !!!subtype_key1),
                       tumor_specific = one_sample$DE_in_pseudobulk,
                       stringsAsFactors = F)
left_anno <- rowAnnotation(df = gen_anno,
                           col = list(subtype_specific= c(subtype_genotype_cols[2:8],
                                                     "Multiple" = "lightgrey"),
                                      tumor_specific = c("Yes" = "Black",
                                                        "No" ="White")))

hm1 <- Heatmap(bulk_norm_expression_mat,
               cluster_columns = T,
               cluster_rows = T,
               show_row_names = F,
               show_column_names  = F,
               top_annotation = top_ano,
               left_annotation = left_anno,
               column_split = anno$Subtype,
               column_title = NULL,
               cluster_row_slices = F,
               cluster_column_slices = F,
               row_title_rot = 0,
               heatmap_legend_param = list(title = "Scaled exp", legend_width = unit(3, "cm"), title_position = "topcenter", direction = "horizontal"))

# dev.off()
# pdf("Figures/bulk_de_heatmap_all_degs.pdf", width =15, height=7)
# hm1
# dev.off()

# ----
# Andrew's pseudobulk z-score expression matrix  
# ----

# Make a heatmap out of the pseudobulk data for normal cells to put next to the bulk
pseudobulk_cpm <- read_csv("/data/gpfs/projects/punim0010/projects/Pattison_projects/PPGL_sc_RNA_seq_compendium/Tables/Pseudobulk_log2_CPM.csv")
pseudobulk_matrix <- as.matrix(pseudobulk_cpm[,2:ncol(pseudobulk_cpm)])
rownames(pseudobulk_matrix) <- pseudobulk_cpm$gene

# Order of cell types that I want
cell_types <- c("Chromaffin.cells", "Sustentacular.cells","Adrenocortical.cells",  "Endothelial.cells", "Fibroblasts", "Myeloid.cells", "T.NK.cells", "B.cells")

# Read in and join on the single nuclei metadata
# Read in the single nuclei metadata
sn_md <- read_tsv("Data/pseudobulk_metadata.tsv")
pseudobulk_annotation <- data.frame(Sample_cell_type = colnames(pseudobulk_matrix),check.names = F)%>%
  mutate(Sample = gsub("_.*" ,"" , Sample_cell_type))%>%
  mutate(Cell_type = gsub(".*_" ,"" ,  Sample_cell_type))%>%
  left_join(sn_md)%>%
  mutate(Cluster = gsub("-", "_", Subgroup3))%>%
  left_join(names_df)%>%
  rename(Subtype = Comete)%>%
  mutate(Cell_type = factor(Cell_type, levels = cell_types))%>%
  arrange(Cell_type)

pseudobulk_matrix <- pseudobulk_matrix[,pseudobulk_annotation$Sample_cell_type]
# Z score transform
# Average over cell types
means_list <- list()
for(i in 1:length(unique(pseudobulk_annotation$Cell_type))){
  cell_type <- unique(pseudobulk_annotation$Cell_type)[i]
  samples <- filter(pseudobulk_annotation, Cell_type == cell_type)
  means <- rowMeans(pseudobulk_matrix[,samples$Sample_cell_type]) %>%
    data.frame(check.rows = F, check.names = F)
  colnames(means) <- cell_type
  means_list[[i]]<- means
}

pseudobulk_matrix <- bind_cols(means_list)%>%
  as.matrix()

pseudobulk_matrix <- t(scale(t(pseudobulk_matrix)))
# Match the rows of the second matrix to the first
matched <- match(rownames(bulk_norm_expression_mat),rownames(pseudobulk_matrix))
pseudobulk_matrix <- pseudobulk_matrix[matched,]
# Set NAs to 0
pseudobulk_matrix[is.na(pseudobulk_matrix)] <- 0
# Reset rownames
rownames(pseudobulk_matrix) == rownames(bulk_norm_expression_mat)

# reformat the column names 
colnames(pseudobulk_matrix) <- gsub(pattern = "\\.",
                                    replacement = " ",
                                    x = colnames(pseudobulk_matrix))
# reformat the column names 
colnames(pseudobulk_matrix) <- gsub(pattern = "T NK",
                                    replacement = "T/NK",
                                    x = colnames(pseudobulk_matrix))

hm2_column_order <- all.cell.types[all.cell.types %in% colnames(pseudobulk_matrix)]

hm2 <- Heatmap(pseudobulk_matrix[,hm2_column_order],
               cluster_columns = FALSE,
               cluster_rows = FALSE,
               column_names_side = "top",
               show_row_names = F,
               show_column_names  = T,
               column_title_rot = 90,
               width = ncol(pseudobulk_matrix[,hm2_column_order])*unit(4, "mm"),
               heatmap_legend_param  = list(title = "Scaled exp", legend_width = unit(3, "cm"), title_position = "topcenter", direction = "horizontal"))

# ----
# Make a single cell heatmap showing fraction expressing for the single nuclei 
# ---- 

# read in the snRNA-seq data 
rna <- readRDS("Data/sn.PPGLs.filtered.with.decontX.RDS")

# add the updated sample metadata 
pb_metadata <- read_tsv("Data/pseudobulk_metadata.tsv") %>% 
  mutate(Genotype_Subgroup4 = paste(Genotype, Subgroup4, sep = "_")) # make a metadata column to aggregate samples with same genotype and subgroup

# check for duplicated column names
dup.cols <- names(pb_metadata)[names(pb_metadata) %in% names(rna[[]])]
# merge the metadata tables
metadata_new <- rna[[]] %>%
  dplyr::select(!dup.cols[2:6]) %>% # remove the old annotations (except sample names) and replace with new ones
  left_join(pb_metadata, by = "orig.ident")
rownames(metadata_new) <- rownames(rna[[]])
rna@meta.data <- metadata_new

# ----
# Calculate fraction cells expressing each gene in each cluster 
# ----

# these were done in two separate functions as there wasn't enough memory to do it in one go
# calculate fraction of normal cells (tumour microenvironment and normal sample)
pseudobulk.fraction.expressed.normal <-  do.call('cbind', lapply(unique(rna$Cell_Type)[unique(rna$Cell_Type) != "Tumour"], function(cell_type){
  rowMeans(as.matrix(rna@assays$RNA@counts[, rna$Cell_Type==cell_type]) > 0)}))
colnames(pseudobulk.fraction.expressed.normal) <- paste(make.names(unique(rna$Cell_Type)[unique(rna$Cell_Type) != "Tumour"]), "fraction.cells.expressing", sep = ".")
# make it a tibble
pseudobulk.fraction.expressed.normal.table <- pseudobulk.fraction.expressed.normal %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  tibble()

# calculate fraction of tumour cells expressing each gene
pseudobulk.fraction.expressed.tumour <-  do.call('cbind', lapply(unique(rna$Subgroup4)[unique(rna$Subgroup4) != "Normal"], function(subgroup){
  rowMeans(as.matrix(rna@assays$RNA@counts[, rna$Subgroup4==subgroup]) > 0)})) # get the number of cells with >0 reads in each row
colnames(pseudobulk.fraction.expressed.tumour) <- paste(make.names(unique(rna$Subgroup4)[unique(rna$Subgroup4) != "Normal"]), "fraction.cells.expressing", sep = ".")
# make it a tibble
pseudobulk.fraction.expressed.tumour.table <- pseudobulk.fraction.expressed.tumour %>%
  data.frame() %>%
  rownames_to_column("gene") %>%
  tibble()

# join the two dataframes
pseudobulk.fraction.expressed.all <- pseudobulk.fraction.expressed.normal.table %>%
  left_join(pseudobulk.fraction.expressed.tumour.table)

#write_tsv(pseudobulk.fraction.expressed.all, "Data/pseudobulk.fraction.expressed.all.tsv")

# ----
# make matrix for the heatmap 
# ---- 

# pseudobulk.fraction.expressed.all <- read_tsv("Data/pseudobulk.fraction.expressed.all.tsv")
# Match the rows of the second matrix to the first, filter out genes I don't want to plot
matched <- match(rownames(bulk_norm_expression_mat),pseudobulk.fraction.expressed.all$gene)
pseudobulk_fraction_cells_expressing <- pseudobulk.fraction.expressed.all[matched,]
# get the gene names
genes_order <- pseudobulk_fraction_cells_expressing$gene
# convert to matrix, remove gene column 
pseudobulk_fraction_cells_expressing_matrix <- as.matrix(pseudobulk_fraction_cells_expressing[2:ncol(pseudobulk_fraction_cells_expressing)])
# set NA counts to zero
pseudobulk_fraction_cells_expressing_matrix[is.na(pseudobulk_fraction_cells_expressing_matrix)] <- 0
# set matrix rownames
rownames(pseudobulk_fraction_cells_expressing_matrix) <- genes_order

# fix the column names 
colnames(pseudobulk_fraction_cells_expressing_matrix) <- gsub(pattern = "\\.fraction\\.cells\\.expressing",
                                         replacement = "",
                                         colnames(pseudobulk_fraction_cells_expressing_matrix))
colnames(pseudobulk_fraction_cells_expressing_matrix) <- gsub(pattern = "\\.",
                                                              replacement = " ",
                                                              colnames(pseudobulk_fraction_cells_expressing_matrix))
# Add the backslash
colnames(pseudobulk_fraction_cells_expressing_matrix)[colnames(pseudobulk_fraction_cells_expressing_matrix) == "Macrophages Monocytes"] <- "Myeloid cells"
colnames(pseudobulk_fraction_cells_expressing_matrix)[colnames(pseudobulk_fraction_cells_expressing_matrix) == "T cells"] <- "T/NK cells"

# specify order of the heatmap columns  
columns.order <- c(
  all.cell.types #,  subgroups.lvl
  )
columns.order <- columns.order[columns.order %in% colnames(pseudobulk_fraction_cells_expressing_matrix)]
# change order of heatmap columns to the right order
pseudobulk_fraction_cells_expressing_matrix <- pseudobulk_fraction_cells_expressing_matrix[,columns.order]

# specify the colour gradient to use for this heatmap
col_fun2 = colorRamp2(
  seq(0.01, 1, len = 11),
  rev(brewer.pal(11, "RdYlBu")))

columns.order
colnames(pseudobulk_fraction_cells_expressing_matrix)

# ----
# organise metadata for annotating the pseudobulk fraction-expressing heatmap
# ----

# make a vector describing if the group is tumour or normal cells to split the heatmap
pb_plot_metadata <- tibble(normal_or_tumour = c(rep("Normal", ncol(pseudobulk.fraction.expressed.normal)),
  rep("Tumour", ncol(pseudobulk.fraction.expressed.tumour))))

# calculate the fraction of total DEGs that are expressed in the normal tissue
# number of DEGs that are expressed in over 20% of cells in each cell type  
ndegs_expresed_in_normal <- colSums(ifelse(pseudobulk_fraction_cells_expressing_matrix > 0.2, 1, 0))
frac_degs_expresed_in_normal <- (ndegs_expresed_in_normal/nrow(pseudobulk_fraction_cells_expressing_matrix))

# ----
# Plot the heatmap 
# ----

# removed mast cells column because no mast cells are present in the pseudobulk expression matrix
hm_deg_fractionexpressing <-  Heatmap(
  pseudobulk_fraction_cells_expressing_matrix[, colnames(pseudobulk_fraction_cells_expressing_matrix)[colnames(pseudobulk_fraction_cells_expressing_matrix) != "Mast cells"]],
  width = ncol(pseudobulk_fraction_cells_expressing_matrix[, colnames(pseudobulk_fraction_cells_expressing_matrix)[colnames(pseudobulk_fraction_cells_expressing_matrix) != "Mast cells"]])*unit(4, "mm"), 
  col = col_fun2,
  column_names_side = "top",
  cluster_rows = FALSE,
  show_row_names = FALSE,
  cluster_columns = FALSE,
  heatmap_legend_param  = list(title = "Cells expressing (%)",
                               legend_width = unit(3, "cm"),
                               title_position = "topcenter",
                               direction = "horizontal"))

hm1 + hm2 + hm_deg_fractionexpressing

dev.off()
pdf("Figures/bulk_de_heatmap_all_degs1.pdf", width = 15, height = 7)
hm1 + hm2 + hm_deg_fractionexpressing
dev.off()