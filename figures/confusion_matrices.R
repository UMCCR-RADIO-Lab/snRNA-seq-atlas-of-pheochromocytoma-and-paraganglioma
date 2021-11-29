# calculating overlap between our classifications
# vs previous subtyping of the same samples
# and visualisation of these results with confusion-matrix-like heatmaps

library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
source("singlecell_colour_palettes.R")

bulk_metadata <- read_csv("Data/Table S3. Pheo-atlas metadata.csv")

bulk_metadata <- bulk_metadata %>%
  filter(Subtype != "Normal") %>% 
  mutate(
    Subgroup4 = recode(Subtype, !!!subtype_key2))

# ----
# make a confusion matrix with the TCGA samples
# ----

tcga_clusters <- bulk_metadata %>%
  select(Sample, Subgroup4, TCGA_Cluster) %>% 
  filter(is.na(TCGA_Cluster) == FALSE) 

# make a matrix showing the intersection of each classification
tcga_cm <- tcga_clusters %>% 
  group_by(Subgroup4, TCGA_Cluster) %>% 
  summarise(n()) %>% 
  ungroup() %>%
  pivot_wider(names_from = TCGA_Cluster, values_from = `n()`) %>% 
  select(Subgroup4, Pseudohypoxia, `Kinase signaling`, `Wnt-altered`, `Cortical admixture`) %>% 
  data.frame()

tcga_cm[is.na(tcga_cm)] <- 0
rownames(tcga_cm) <- tcga_cm$Subgroup4
tcga_cm$Subgroup4 <- NULL
tcga_cm <- as.matrix(tcga_cm)

col_fun = colorRamp2(
  seq(0, 80, len = 9),
  brewer.pal(9, "Blues"))

tcga_hm <- Heatmap(tcga_cm,
        cluster_columns = F,
        column_names_side = "top",
        col = col_fun,
        cluster_rows = F, 
        row_names_side = "left",
        width = ncol(tcga_cm)*unit(20, "mm"), 
        height = nrow(tcga_cm)*unit(20, "mm"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.0f", tcga_cm[i, j]), x, y, gp = gpar(fontsize = 10))},
        #heatmap_legend_param  = list(title = "")
        show_heatmap_legend = F)

# ----
# Castro-vega clusters confusion matrix
# ----

castrovega_clusters <- bulk_metadata %>%
  select(Sample, Subgroup4, Castro_Vega_mRNA_Classification)

# make a matrix showing the intersection of each classification
castrovega_cm <- castrovega_clusters %>% 
  group_by(Subgroup4, Castro_Vega_mRNA_Classification) %>% 
  summarise(n()) %>% 
  ungroup() %>%
  pivot_wider(names_from = Castro_Vega_mRNA_Classification, values_from = `n()`) %>% 
  select(Subgroup4, C1A, C1B, C2A, C2B, C2C) %>% 
  data.frame()

castrovega_cm[is.na(castrovega_cm)] <- 0
rownames(castrovega_cm) <- castrovega_cm$Subgroup4
castrovega_cm$Subgroup4 <- NULL
castrovega_cm <- as.matrix(castrovega_cm)

col_fun = colorRamp2(
  seq(0, 80, len = 9),
  brewer.pal(9, "Blues"))

castrovega_hm <- Heatmap(castrovega_cm,
                   cluster_columns = F,
                   column_names_side = "top",
                   col = col_fun,
                   cluster_rows = F, 
                   row_names_side = "left",
                   width = ncol(castrovega_cm)*unit(20, "mm"), 
                   height = nrow(castrovega_cm)*unit(20, "mm"), 
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%.0f", castrovega_cm[i, j]), x, y, gp = gpar(fontsize = 10))},
                   #heatmap_legend_param  = list(title = "")
                   show_heatmap_legend = F)

castrovega_hm + tcga_hm

# dev.off()
# pdf(file = "Figures/confusion_matrices.pdf", height = 10, width = 10)
# castrovega_hm + tcga_hm
# dev.off()


