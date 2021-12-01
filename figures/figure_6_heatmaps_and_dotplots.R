rm(list=ls())

library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(patchwork)
library(ggsci)
library(scales)
library(ComplexHeatmap)
library(circlize)
library(org.Hs.eg.db)
source("Figures/dotplot_functions.R")
source("singlecell_colour_palettes.R")

# ----
# Organise Jansky data
# ----

# read in the Jansky data
fetal_medulla_jansky <- readRDS("Data/Jansky_2021/snRNA_adrenal_medulla_jansky.RDS")
fetal_medulla_jansky$cell_type <- Idents(object = fetal_medulla_jansky)

discrete_palette <- pal_d3(palette = "category20")(length(unique(Idents(fetal_medulla_jansky))))
colours_jansky <- setNames(object = discrete_palette, nm = as.character(unique(Idents(fetal_medulla_jansky))))

# remove the SCP cell types
fetal_medulla_jansky <- subset(fetal_medulla_jansky,
                               subset = cell_type %in% c("SCPs","cycling SCPs", "late SCPs" ), 
                               invert = TRUE)

jansky_md <- fetal_medulla_jansky@meta.data

jansky_md <- jansky_md %>%
  mutate(cell_type_plot = recode(cell_type,
                                 "connecting Chromaffin cells" = "Connecting progenitor cells",
                                 "Chromaffin cells" = "Early chromaffin cells",
                                 "late Chromaffin cells" = "Late chromaffin cells",
                                 "Neuroblasts" = "Early neuroblasts",
                                 "cycling Neuroblasts" = "Cycling neuroblasts",
                                 "late Neuroblasts" = "Late neuroblasts"))
fetal_medulla_jansky@meta.data <- jansky_md

# ----
# Organise the single cell PCPG data 
# ----

# read in the PCPG data
plot_rna <- readRDS("Data/pcpg_chromaffincells.RDS")

# remove the SCLCs
plot_rna <- subset(plot_rna,
                   subset = Subgroup4 %in% c("Tumour-associated SLCs","Normal SLCs"),
                   invert=TRUE)

pcpg_md <- plot_rna@meta.data
subtypes_genotypes <- c("CCs (NAM)", "C1A1 (SDHx)", "C1A2 (SDHx-HN)", "C1B1 (VHL)", "C1B2 (EPAS1)", "C2A (Kinase)", "C2B1 (MAX)", "C2B2 (MAML3)", "C2C")
subtypes_old <- c("Chromaffin cells", "C1Ai", "C1Aii", "C1Bi", "C1Bii", "C2A", "C2Bi", "C2Bii", "C2C")
subtype_key1 = setNames(subtypes_genotypes, subtypes_old)
pcpg_md <- pcpg_md %>%
  # Recode the subtype for plotting
  mutate(Subgroup4 = recode(Subgroup4,!!!subtype_key1)) %>% 
  mutate(Genotype = if_else(Subgroup4 == "CCs (NAM)","Normal", Genotype))
# Recode the genotype for plotting

plot_rna@meta.data <- pcpg_md

# manually arrange genotype order
plot_rna$Genotype <-  factor(plot_rna$Genotype, levels = c("SDHA", "SDHB", "SDHD", "VHL", "EPAS1",  "HRAS", "NF1", "RET", "TMEM127", "H3F3A", "MAX", "Unknown", "FH", "MAML3", "Normal"))

# ----
# Make a list of  fetal marker genes from Jansky, and additional interesting genes 
# ----

# TODO: mark the receptors in adobe illustrator

# medulla markers
bridge_cell_markers <- c("ASCL1",
                         "CDH9",
                         "ERBB4")

chromaffin_markers <- c("TH",
                        "PENK",
                        "CARTPT",
                        "MEG3",
                        "DLK1",
                        "PNMT")

neuroblast_markers <- c("ALK",
                        "PHOX2A",
                        "NTRK3",
                        "RET",
                        "NPY")
hypoxia <- c("EPAS1",
             "VEGFA",
             "EGLN1",
             "EGLN3",
             "SLC2A1")

C1Bi_up <- c("SSTR2")

C1Bi_up <- c("GPR139",
             "PTHLH",
             "EGFR",
             "TFAP2C",
             "AQP1")

C2Bii_up <- c("IRX4",
              "HMGA2",
              "HMX1",
              "AQP2",
              "SST",
              "SSTR2",
              "AR",
              "KIT",
              "GLI2",
              "WNT4")

MAML3_fusion <- c("MAML3",
                  "UBTF",
                  "TCF4")

genes_plot <- c(
  bridge_cell_markers,
  chromaffin_markers,
  neuroblast_markers,
  hypoxia,
  C1Bi_up, 
  C2Bii_up,
  MAML3_fusion)

gene_ontologies <- c(
  rep("Bridge cell", length(bridge_cell_markers)),
  rep("Chromaffin cell", length(chromaffin_markers)),
  rep("Neuroblast", length(neuroblast_markers)),
  rep("Hypoxia", length(hypoxia)),
  rep("C1B1 (VHL) up", length(C1Bi_up)),
  rep("C2B2 (MAML3) up", length(C2Bii_up)),
  rep("MAML3 fusion", length(MAML3_fusion)))

gene_ontologies <- factor(gene_ontologies,
                          levels = c("Bridge cell",
                                     "Chromaffin cell",
                                     "Neuroblast",
                                     "Hypoxia",
                                     "C1B1 (VHL) up",
                                     "C2B2 (MAML3) up",
                                     "MAML3 fusion"))

genes_plot_df <- tibble(gene = genes_plot,
                        gene_ontology = gene_ontologies)

# ----
# Organise the bulk data 
# ---- 

# specify plotting order
genotypes.lvl <- c("Normal", "SDHB", "SDHA", "SDHD", "VHL", "EPAS1", "FH", "HRAS", "NF1", "RET", "TMEM127", "MAX", "H3F3A", "Unknown", "MAML3")

# read in the bulk rna-seq metadata
bulk_metadata <- read_tsv("Data/bulk_metadata_gsva.tsv")
# specify order I want the subgroups to be plotted
bulk_metadata <- bulk_metadata %>%
  mutate(Subgroup4 = recode(Subgroup4, !!!subtype_key1))
bulk_metadata$Subgroup4 <- factor(bulk_metadata$Subgroup4, levels = c(subtypes_genotypes))

# Read in bulk data (batch-normalised expression matrix)
bulk_rna <- read_tsv("Data/bulk_expression_batch_normalised.tsv")

# remove the samples from the bulk that aren't in the metadata (normals and adrenocortical admixture)
bulk_rna <- bulk_rna %>%
  dplyr::select(all_of(c("Gene", bulk_metadata$Sample_raw)))

# Get the batch normed expression from the bulk
bulk_plot_mat <- as.matrix(bulk_rna[,2:ncol(bulk_rna)])
# Add genes as rownames
rownames(bulk_plot_mat) <- bulk_rna$Gene
# Z score transform
bulk_plot_mat <- t(scale(t(bulk_plot_mat)))

# filter out the genes that arent in the bulk (1 LncRNA gene)
genes_plot_df <- genes_plot_df[genes_plot_df$gene %in% rownames(bulk_plot_mat), ]
# subset the gene expression matrix to just include the genes i want to plot
bulk_heatmap_deg <- bulk_plot_mat[genes_plot_df$gene, ]

# organise heatmap annotations 
# check the order of samples in the metadata matches the gex matrix 
table(bulk_metadata$Sample == colnames(bulk_plot_mat))


# ----
# Plot curated gene expression in the pcpg
# ----

dot_deg_pcpg = HeatmapDotPlot.Seurat(plot_rna,
                                     #split.by = "Genotype", 
                                     features = genes_plot_df$gene,
                                     gene_grouping = genes_plot_df$gene_ontology,
                                     aggr.by = c("Genotype","Subgroup4"),
                                     split.by = "Subgroup4",
                                     column_title = NULL,
                                     annot.columns = c("Subgroup4"),
                                     annot.colours = list("Subgroup4" = subtype_genotype_cols),
                                     assay = "SCT",
                                     show_annotation_name = FALSE,
                                     slot = "data",
                                     col = c("blue", "grey", "red"),
                                     column_title_rot = 90,
                                     row_names_gp = gpar(fontface = "italic"),
                                     # left_annotation = side_anno_deg,
                                     show_column_names = TRUE,
                                     column_names_side = "top",
                                     cluster_columns = FALSE,
                                     cluster_column_slices = FALSE,
                                     cluster_rows = FALSE,
                                     show_column_dend = FALSE,
                                     show_row_dend = TRUE,
                                     show_legend = FALSE,
                                     heatmap_legend_param  = hm_legend_params )

# dev.off()
# pdf(file = "Figures/pcpg_dotplot_degs.pdf", height = 28, width = 18)
# dot_deg_pcpg
# dev.off()

# double check that the samples and annotation are in the same order 
#colnames(dot_deg_pcpg@matrix) == sample.order

# ----
# Plot the curated genes in the Jansky Data 
# ---- 

dot_deg_fetal = HeatmapDotPlot.Seurat(fetal_medulla_jansky,
                                      features = genes_plot_df$gene,
                                      gene_grouping = genes_plot_df$gene_ontology,
                                      aggr.by = "cell_type_plot",
                                      cluster_row_slices = FALSE,
                                      cluster_columns = FALSE,
                                      cluster_rows = FALSE,
                                      assay = "RNA",
                                      row_names_side = "left",
                                      column_names_side = "top",
                                      slot = "data",
                                      col = c("blue", "grey", "red"),
                                      column_title_rot = 90,
                                      # left_annotation = side_anno_deg,
                                      row_names_gp = gpar(fontface = "italic"),
                                      heatmap_legend_param  = hm_legend_params)

# ----
# plot the curated genes in the Bulk
# ----

# make an annotation with all the GSVA scores 
# to get the 1 and 99th percentile cat all the gsva scores cat these all together 
gsva_all_scores <- c(bulk_metadata$gsva_bridge_cell,
                     bulk_metadata$gsva_connecting_chromaffin_cell,
                     bulk_metadata$gsva_chromaffin_cell,
                     bulk_metadata$gsva_late_chromaffin_cell,
                     bulk_metadata$gsva_neuroblast,
                     bulk_metadata$gsva_late_neuroblast)
# this function automatically sets the scale between the 0.01th and 99th percentile
greys_col_fun = colorRamp2(seq(quantile(gsva_all_scores, 0.01, na.rm = TRUE), quantile(gsva_all_scores, 0.99, na.rm = TRUE), len = 9), brewer.pal(9, "Greys"))

top_anno_bulk <- HeatmapAnnotation(
  "Subtype" = bulk_metadata$Subgroup4,
  #"Genotype" = bulk_metadata$Genotype,
  "Bridge cell GSVA" = bulk_metadata$gsva_bridge_cell,
  "Connecting progenitor GSVA" = bulk_metadata$gsva_connecting_chromaffin_cell,
  "Early chromaffin cell GSVA" = bulk_metadata$gsva_chromaffin_cell,
  "Late chromaffin cell GSVA" = bulk_metadata$gsva_late_chromaffin_cell,
  "Neuroblast GSVA" = bulk_metadata$gsva_neuroblast,
  "Late neuroblast GSVA" = bulk_metadata$gsva_late_neuroblast,
  show_annotation_name = c("Subtype" = FALSE),
  annotation_legend_param = list(
    "Subtype" = list(legend_height = unit(3, "cm")),
    # use this legend for all gsva scores
    "Bridge cell GSVA" = list(title = "GSVA score",
                              # title_position = "topleft",
                              legend_width = unit(3, "cm"),
                              legend_height = unit(3, "cm")
                              # direction = "vertical",
                              # border = "grey"
    )),
  col = list("Subtype" = subtype_genotype_cols,
             #"Genotype" = genotype.cols,
             "Bridge cell GSVA" = greys_col_fun,
             "Connecting progenitor GSVA" = greys_col_fun,
             "Early chromaffin cell GSVA" = greys_col_fun,
             "Late chromaffin cell GSVA" = greys_col_fun,
             "Neuroblast GSVA" = greys_col_fun, 
             "Late neuroblast GSVA" = greys_col_fun), 
  show_legend = c("Subtype" = TRUE,
                  "Genotype" = FALSE,
                  "Bridge cell GSVA" = TRUE, 
                  "Connecting progenitor GSVA" = FALSE, 
                  "Early chromaffin cell GSVA" = FALSE, 
                  "Late chromaffin cell GSVA" = FALSE, 
                  "Neuroblast GSVA" = FALSE, 
                  "Late neuroblast GSVA" = FALSE))

# plot the heatmap 
bulk_hm_deg <- Heatmap(bulk_heatmap_deg,
                       # row_split = genes.plot$Cluster,
                       column_split = bulk_metadata$Subgroup4,
                       width = ncol(bulk_heatmap_deg)*unit(0.3, "mm"), 
                       height = unit(300, "mm"),
                       row_gap = unit(2, "mm"), 
                       column_gap = unit(2, "mm"),
                       row_title_rot = 0,
                       top_annotation = top_anno_bulk,
                       column_title_rot = 90,
                       cluster_columns = TRUE, 
                       show_column_dend = FALSE,
                       cluster_column_slices = FALSE,
                       cluster_rows = FALSE,
                       show_row_dend = FALSE,
                       cluster_row_slices = FALSE,
                       show_row_names = TRUE,
                       row_names_side = "right",
                       show_column_names  = FALSE,
                       row_names_gp = gpar(fontface = "italic"),
                       heatmap_legend_param  = hm_legend_params)

cell.size = 0.5
label="Scaled exp"
size.label="Fraction expressing"
hm.legend = list()
hm.legend = c(hm.legend, list(Legend(title = size.label,
                                     labels = c(0.25, 0.50, 0.75, 1.00) %>% as.character,
                                     size=unit.c(unit(sqrt(0.25)*cell.size, "cm"),
                                                 unit(sqrt(0.5)*cell.size, "cm"),
                                                 unit(sqrt(0.75)*cell.size , "cm"),
                                                 unit(sqrt(1.0)*cell.size, "cm")),
                                     type = "points",
                                     grid_height = unit(cell.size,"cm"),
                                     grid_width=unit(cell.size,"cm"),
                                     legend_height=unit(4*cell.size*2, "cm"),
                                     background = NA)))

# ----
# Combine the plots into one 
# ---- 

dev.off()
pdf(file = "Figures/jansky_dotplot_degs_genotype.pdf",
    height = 28, width = 24)
draw(dot_deg_fetal + dot_deg_pcpg + bulk_hm_deg, annotation_legend_list=hm.legend,
     merge_legend = TRUE,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "bottom")
dev.off()


# ----
# HIF gene correlation matrix 
# ----

# The HIF target genes are (mostly) pheo-associated HIF signalling genes from Jochmanova et al (https://academic.oup.com/jnci/article/105/17/1270/908345?login=true)
# A few were filtered out from the single cell data (presumably due to low expression in all cell types?)

# this matrix was made in the HIF_target_correlation.R script
spearman_correlations <- read.csv("Results/VEGFA_spearman_correlation_HIF_target_genes.tsv",
                                  row.names = 1)
spearman_correlations_mat <- as.matrix(spearman_correlations)
# these genes were removed as they had low expression and didn't correlate with the other genes
genes.remove <- c("CCND1",
                  "TGFB3",
                  "NOS2",
                  "ALDOA",
                  "ADM",
                  "PFKM")
spearman_correlations_mat <- spearman_correlations_mat[!rownames(spearman_correlations_mat) %in% genes.remove,
                                                       !colnames(spearman_correlations_mat) %in% genes.remove]

ylgnbu_col_fun = colorRamp2(seq(-1,1,len = 9), rev(brewer.pal(9, "RdYlBu")))

# plot the Spearman correlation heatmap 
hm_hif_target_cor_spearman <- Heatmap(spearman_correlations_mat,
                                      width = ncol(spearman_correlations_mat)*unit(5, "mm"), 
                                      height = nrow(spearman_correlations_mat)*unit(5, "mm"),
                                      col = ylgnbu_col_fun,
                                      column_title_rot = 90,
                                      row_names_side = "left",
                                      row_dend_side = "left",
                                      row_names_gp = gpar(fontface = "italic"),
                                      column_names_gp = gpar(fontface = "italic"),
                                      heatmap_legend_param  = list(
                                        title = "Rho", 
                                        legend_width = unit(3, "cm"),
                                        title_position = "topcenter",
                                        direction = "horizontal"))
hm_hif_target_cor_spearman

# ----
# plot HIF target genes across the subtypes 
# ----

HIF_genes <- c(rownames(hm_hif_target_cor_spearman@matrix))

# ----
# HIF genes in FETAL 
# ----

dot_hif_fetal = HeatmapDotPlot.Seurat(fetal_medulla_jansky,
                                      features = HIF_genes,
                                      # gene_grouping = genes_plot_df$gene_ontology,
                                      aggr.by = "cell_type_plot",
                                      cluster_row_slices = FALSE,
                                      cluster_columns = FALSE,
                                      cluster_rows = FALSE,
                                      assay = "RNA",
                                      row_names_side = "left",
                                      column_names_side = "top",
                                      slot = "data",
                                      col = c("blue", "grey", "red"),
                                      column_title_rot = 90,
                                      # left_annotation = side_anno_deg,
                                      row_names_gp = gpar(fontface = "italic"),
                                      heatmap_legend_param  = hm_legend_params)

dot_hif_pcpg = HeatmapDotPlot.Seurat(plot_rna,
                                     split.by = "Subgroup4", 
                                     features = HIF_genes,
                                     aggr.by = c("Genotype", "Subgroup4"),
                                     annot.columns = c("Subgroup4"),
                                     annot.colours = list("Subgroup4" = subtype_genotype_cols),
                                     assay = "SCT",
                                     show_annotation_name = FALSE,
                                     slot = "data",
                                     col = c("blue", "grey", "red"),
                                     column_title_rot = 90,
                                     row_names_gp = gpar(fontface = "italic"),
                                     # left_annotation = side_anno_deg,
                                     show_column_names = TRUE,
                                     column_names_side = 'top',
                                     cluster_columns = FALSE,
                                     cluster_column_slices = FALSE,
                                     cluster_rows = FALSE,
                                     show_column_dend = FALSE,
                                     show_row_dend = TRUE,
                                     show_legend = TRUE,
                                     heatmap_legend_param = hm_legend_params)

# ----
# HIF genes in BULK
# ----

# subset the gene expression matrix to just include the genes i want to plot
bulk_heatmap_hif <- bulk_plot_mat[HIF_genes, ]

# make an annotation wihtout the GSVA part
top_anno_bulk_no_gsva <- HeatmapAnnotation(
  "Subtype" = bulk_metadata$Subgroup4,
  show_annotation_name = c("Subtype" = FALSE),
  annotation_legend_param = list(
    "Subtype" = list(legend_height = unit(3, "cm"))),
  col = list("Subtype" = cluster.cols), 
  show_legend = c("Subtype" = TRUE))

# plot the heatmap 
bulk_hm_hif <- Heatmap(bulk_heatmap_hif,
                       # row_split = genes.plot$Cluster,
                       column_split = bulk_metadata$Subgroup4,
                       width = ncol(bulk_heatmap_deg)*unit(0.3, "mm"), 
                       height = unit(300, "mm"),
                       row_gap = unit(2, "mm"), 
                       column_gap = unit(2, "mm"),
                       row_title_rot = 0,
                       top_annotation = top_anno_bulk_no_gsva,
                       column_title_rot = 90,
                       cluster_columns = TRUE, 
                       show_column_dend = FALSE,
                       cluster_column_slices = FALSE,
                       cluster_rows = FALSE,
                       show_row_dend = FALSE,
                       cluster_row_slices = FALSE,
                       show_row_names = TRUE,
                       row_names_side = "right",
                       show_column_names  = FALSE,
                       row_names_gp = gpar(fontface = "italic"),
                       heatmap_legend_param  = hm_legend_params)

# dev.off()
# pdf(file = "Figures/HIF_genes_dotplot.pdf", height = 24, width = 30)
# draw(hm_hif_target_cor_spearman + dot_hif_pcpg, 
#      annotation_legend_list=hm.legend, 
#      merge_legend = TRUE,
#      heatmap_legend_side = "bottom",
#      annotation_legend_side = "bottom")
# dev.off()

