# Analysis of the PCPG cohort using harmony integration for batch correction

rm(list = ls())

library(Seurat)
library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(scales)
library(speckle)
library(harmony)
library(circlize)
library(Matrix.utils)

source("singlecell_utility_functions.R")
source("Figures/dotplot_functions.R")
source("singlecell_colour_palettes.R")

# ----
# read and prep data
# ----

# Seurat object with single cell data
pcpg_rna <- readRDS("Data/pcpg_with_metadata_and_qc.RDS") 

uncorrected_umaps <- DimPlot(pcpg_rna, group.by = "Patient") +
  DimPlot(pcpg_rna, group.by = "tumor_subtype") +
  DimPlot(pcpg_rna, group.by = "Sample")

# looks like a strong patient-specific effect in the tumor cells - this can be removed with RPCA batch correction
uncorrected_umaps

md <- pcpg_rna@meta.data

# ----
# batch correction with harmony
# ----

DefaultAssay(pcpg_rna) <- "RNA"
pcpg_rna@assays$SCT <- NULL # remove the old SCT assay to be safe

# run the standard Seurat workflow up to PCA step
pcpg_rna_harmonized <- CreateSeuratObject(counts = pcpg_rna@assays$RNA@counts,
                                          meta.data = md,
                                          project = "pcpg_snRNA_seq") %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method="vst", nfeatures=3000, verbose=FALSE) %>% 
  ScaleData() %>% 
  # ScaleData(verbose = FALSE, vars.to.regress=c('S.Score', 'G2M.Score')) %>% 
  RunPCA(verbose = FALSE) 

# Run harmony to remove patient-specific effects
pcpg_rna_harmonized <- RunHarmony(pcpg_rna_harmonized, group.by.vars="Patient", plot_convergence = TRUE)

pcpg_rna_harmonized <- pcpg_rna_harmonized %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5)

saveRDS(pcpg_rna_harmonized, "Results/pcpg_harmonized_patient.RDS")

# ----
# Plot the corrected data 
# ----

# pcpg_rna_harmonized <- readRDS("Results/pcpg_harmonized_patient.RDS")

patient_umap <- DimPlot(pcpg_rna_harmonized, group.by = "Patient") + NoLegend()
batch_umap <- DimPlot(pcpg_rna_harmonized, group.by = "Processing_batch")
cell_types_umap <- DimPlot(pcpg_rna_harmonized, group.by = "Cell_Type")
subtypes_umap <- DimPlot(pcpg_rna_harmonized, group.by = "tumor_subtype")
clusters_umap <-  DimPlot(pcpg_rna_harmonized, group.by =  "seurat_clusters") + NoLegend()
phase_umap <-  DimPlot(pcpg_rna_harmonized, group.by =  "Phase")
types_umap <-  DimPlot(pcpg_rna_harmonized, group.by =  "pseudohypoxic_or_kinase_signalling")
genotypes_umap <- DimPlot(subset(pcpg_rna_harmonized, subset = Cell_Type %in% c("Tumour","Chromaffin cells")), group.by = "Genotype")
phase_umap <- DimPlot(pcpg_rna_harmonized, group.by = "Phase")

subtypes_umap_split <- DimPlot(subset(pcpg_rna_harmonized,
                                      subset = Cell_Type %in% c("Tumour","Chromaffin cells")),
                               group.by = "Sample",
                               split.by="tumor_subtype") + NoLegend()

umaps <- (cell_types_umap + patient_umap + batch_umap +
    types_umap + subtypes_umap + phase_umap) 

ggsave(plot = umaps, filename = "Figures/harmonised_umaps.pdf",
       width = unit(15, "cm"),
       height = unit(7.5, "cm"))

marker_genes <- c("PNMT",
                  "TH",
                  "CHGA",
                  "NTRK3",
                  "LEF1",
                  "EPAS1")

featureplots <- FeaturePlot(pcpg_rna_harmonized,
            features = marker_genes,
            order = TRUE)
featureplots <- featureplots + plot_layout(ncol=3)

ggsave(plot = featureplots, filename = "Figures/harmonised_featureplots.pdf",
       width = unit(15, "cm"),
       height = unit(7.5, "cm"))


# ########################################################################
# Tumour cells subclustering
# ########################################################################


# ----
# TODO:investigate intratumoral heterogeneity
# ----

# get the tumor cells only and reanalyse them 

# MODULE SCORES ---- 

FeaturePlot(pcpg_rna_harmonized,
            features = colnames(module_scores)[1:6],
            order = TRUE)
FeaturePlot(pcpg_rna_harmonized,
            features = colnames(module_scores)[7:12],
            order = TRUE)

# SEURAT CLUSTERS ----

# find markers for all of the tumor clusters
# markers <- FindAllMarkers(pcpg_rna_harmonized)
# write_csv(markers, "Results/tumor_batch_corrected_markers_seurat.csv")
# 
# markers_filtered <- markers %>%
#   filter(p_val_adj < 0.05) %>%
#   group_by(cluster) %>%
#   slice_max(n=5, order_by= avg_logFC)
# 
# FeaturePlot(pcpg_rna_harmonized, features = c("NTNG1", "CARTPT", "CHGA"))

# ########################################################################
# Macrophages subclustering
# ########################################################################

# ----
# subcluster the macrophages
# ----

# get counts for macrophages only
macrophage_rna <- subset(pcpg_rna_harmonized,
                         Cell_Type == "Myeloid cells" &
                           cell_subtype == "Macrophages")
# check that this has just kept all macs
DimPlot(macrophage_rna, group.by = "cell_subtype")

# rerun Seurat pipeline 
macrophage_rna <- macrophage_rna %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method="vst", nfeatures=2000, verbose=FALSE) %>% 
  ScaleData(verbose = FALSE, vars.to.regress=c('S.Score', 'G2M.Score')) %>% 
  RunPCA(verbose = FALSE) 

# print out the PCA Dimensions
ProjectDim(macrophage_rna, reduction = "pca", dims.print = 1:5) 
# chromaffin and adrenocortical cell signature genes and XIST are in top dimensions, indicating ambient RNA and gender-derived clustering

# PCA plots
DimPlot(macrophage_rna, group.by = "Patient", reduction="pca")
# looks like there is clustering by sample in the PCA plot

# run harmony to remove patient-specific effects
macrophage_rna <- RunHarmony(macrophage_rna, group.by.vars="Patient", plot_convergence = TRUE)

# PCA plot of corrected data
DimPlot(macrophage_rna, group.by = "Patient", reduction="harmony")
# no more sample-specific clustering

# downstream analysis
macrophage_rna <- macrophage_rna %>% 
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)

# print the Harmony PCA dimensions
ProjectDim(macrophage_rna, reduction = "harmony", dims.print = 1:5) 
# Explore heatmap of PCs
DimHeatmap(macrophage_rna, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

# look for batch / patient / cell cycle derived clustering
(DimPlot(macrophage_rna, group.by = "Patient") + NoLegend()) +
  DimPlot(macrophage_rna, group.by = "Processing_batch")+ 
  # look at biological factors
  DimPlot(macrophage_rna, group.by = "Phase") +
  (DimPlot(macrophage_rna, group.by = "seurat_clusters") + NoLegend()) +
  DimPlot(macrophage_rna, group.by = "tumor_subtype") +
  DimPlot(macrophage_rna, group.by = "Genotype") +
  DimPlot(macrophage_rna, group.by = "tumor_or_normal")

# ----
# look at the cluster abundance for each tumor subtype 
# ----

# find markers for all of the macrophage clusters
mac_markers <- FindAllMarkers(macrophage_rna)
write_csv(mac_markers, "Results/macrophage_cluster_markers_seurat.csv")
cluster_umap <- DimPlot(macrophage_rna, group.by = 'seurat_clusters',
                        cols = DiscretePalette(n = length(
                          unique(macrophage_rna$seurat_clusters))))

mac_markers <- read_csv("Results/macrophage_cluster_markers_seurat.csv")
markers_filtered <- mac_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by= avg_logFC)
# heatmap of the top 5 marker genes of each cluster
cluster_heatmap <- DoHeatmap(macrophage_rna, features=markers_filtered$gene,
          group.colors =  DiscretePalette(n = length(
            unique(macrophage_rna$seurat_clusters))))

cluster_umap + cluster_heatmap + plot_annotation(title = 'Macrophages')

pal <- setNames(DiscretePalette(n = 12), nm = as.character(0:11))
# plot the abundances
# cluster abundance proportions of samples 
mac_clusters_proportions <- macrophage_rna@meta.data %>%
  dplyr::select(seurat_clusters,
                Sample,
                Genotype,
                tumor_subtype) %>%
  group_by(seurat_clusters,
           Sample,
           Genotype,
           tumor_subtype) %>%
  summarise(cluster_count = n()) %>% # cell count in groups
  group_by(Sample) %>% 
  mutate(cluster_total = sum(cluster_count)) %>%
  mutate(cluster_prop = cluster_count / cluster_total)

mac_clusters_proportions <- ggplot(mac_clusters_proportions,
                                   aes(x = Sample,
                                       y = cluster_prop,
                                       fill = seurat_clusters)) +
  geom_col(position = 'fill') +
  theme_classic() +
  facet_grid(~tumor_subtype + Genotype, scales = "free_x", space = "free_x")+
  labs(y = "Fraction") +
  scale_fill_manual(values =  pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank())

mac_clusters_proportions +
  labs(colour = "") +
  plot_layout(guides = 'collect') 

# featureplots
FeaturePlot(macrophage_rna, features = c("TNF", "EGR2", "IL1B", "EGR1", "HLA-DRA", "CD86"), order = TRUE)
# 
inflammatory <- c("TNF", "TNFSF18", "EGR2", "IL1B", "CCL4", "CCL3", "CCL2", "PLAU", "PLAUR")
M1 <- c("HLA-DRA", "CD86", "CD80")
M2 <- c("CD163", "CD68", "MRC1")
FeaturePlot(macrophage_rna, features = c(inflammatory, M1, M2), order = TRUE)

# ----
# pseudobulk dge for the macrophage subclusters
# ----
# 
# macrophage_cluster_pseudobulk <- make_pseudobulk_matrix(obj = macrophage_rna, grp = "seurat_clusters")
# 
# markers_filtered2 <- mac_markers %>%
#   filter(p_val_adj < 0.05) %>%
#   filter(cluster %in% c(3,9)) %>% 
#   group_by(cluster) %>%
#   slice_max(n=20, order_by= avg_logFC)
# cluster_heatmap2 <- DoHeatmap(macrophage_rna, features=markers_filtered2$gene,
#                              group.colors =  DiscretePalette(n = length(
#                                unique(macrophage_rna$seurat_clusters))))


# ########################################################################
# Angiogenic subclustering
# ########################################################################


# -----
# FIBROBLASTS
# -----

# pcpg_rna_harmonized <- readRDS("Results/pcpg_harmonized_patient.RDS")

# get counts for angiogenic cell types only
fibroblasts_rna <- subset(pcpg_rna_harmonized,
                         Cell_Type == "Fibroblasts")
# check that this has just kept all fibs
DimPlot(fibroblasts_rna, group.by = "cell_subtype")

# rerun Seurat pipeline 
fibroblasts_rna <- fibroblasts_rna %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method="vst", nfeatures=2000, verbose=FALSE) %>% 
  ScaleData(verbose = TRUE, vars.to.regress=c('S.Score', 'G2M.Score')) %>% 
  RunPCA(verbose = FALSE) 

# print out the PCA Dimensions
ProjectDim(fibroblasts_rna, reduction = "pca", dims.print = 1:5) 
# chromaffin and adrenocortical cell signature genes and XIST are in top dimensions, indicating ambient RNA and gender-derived clustering

# PCA plots
DimPlot(fibroblasts_rna, group.by = "Patient", reduction="pca")
# looks like there is clustering by sample in the PCA plot

# run harmony to remove patient-specific effects
fibroblasts_rna <- RunHarmony(fibroblasts_rna, group.by.vars="Patient", plot_convergence = TRUE)

# PCA plot of corrected data
DimPlot(fibroblasts_rna, group.by = "Patient", reduction="harmony")
# no more sample-specific clustering

# downstream analysis
fibroblasts_rna <- fibroblasts_rna %>% 
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)

# print the Harmony PCA dimensions
ProjectDim(fibroblasts_rna, reduction = "harmony", dims.print = 1:5) 
# Explore heatmap of PCs
DimHeatmap(fibroblasts_rna, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

DimPlot(fibroblasts_rna, group.by = 'cell_subtype') +
  DimPlot(fibroblasts_rna, group.by = 'tumor_or_normal') + 
  DimPlot(fibroblasts_rna, group.by = 'Patient') +
  DimPlot(fibroblasts_rna, group.by = 'Phase') +
  DimPlot(fibroblasts_rna, group.by = 'seurat_clusters')

fib_markers <- FindAllMarkers(fibroblasts_rna)
write_csv(fib_markers, "Results/fibroblasts_cluster_markers_seurat.csv")

fibroblasts_cluster_umap <- DimPlot(fibroblasts_rna, group.by = 'seurat_clusters',
                                    cols = DiscretePalette(n = length(
                                      unique(fibroblasts_rna$seurat_clusters))))

fib_markers <- read_csv("Results/fibroblasts_cluster_markers_seurat.csv")
fibroblasts_markers_filtered <- fib_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by= avg_logFC)
# heatmap of the top 5 marker genes of each cluster
fibroblasts_cluster_heatmap <- DoHeatmap(fibroblasts_rna, features=fibroblasts_markers_filtered$gene,
                                         group.colors =  DiscretePalette(n = length(
                                           unique(fibroblasts_rna$seurat_clusters))))

fibroblasts_cluster_umap + fibroblasts_cluster_heatmap + plot_annotation(title = 'fibroblasts')

# -----
# ENDOTHELIAL CELLS
# -----

# get counts for angiogenic cell types only
endothelial_rna <- subset(pcpg_rna_harmonized,
                          Cell_Type == "Endothelial cells")
# check that this has just kept all fibs
DimPlot(endothelial_rna, group.by = "cell_subtype")

# rerun Seurat pipeline 
endothelial_rna <- endothelial_rna %>% 
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method="vst", nfeatures=2000, verbose=FALSE) %>% 
  ScaleData(verbose = TRUE, vars.to.regress=c('S.Score', 'G2M.Score')) %>% 
  RunPCA(verbose = FALSE) 

# print out the PCA Dimensions
ProjectDim(endothelial_rna, reduction = "pca", dims.print = 1:5) 
# chromaffin and adrenocortical cell signature genes and XIST are in top dimensions, indicating ambient RNA and gender-derived clustering

# PCA plots
DimPlot(endothelial_rna, group.by = "Patient", reduction="pca")
# looks like there is clustering by sample in the PCA plot

# run harmony to remove patient-specific effects
endothelial_rna <- RunHarmony(endothelial_rna, group.by.vars="Patient", plot_convergence = TRUE)

# PCA plot of corrected data
DimPlot(endothelial_rna, group.by = "Patient", reduction="harmony")
# no more sample-specific clustering

# downstream analysis
endothelial_rna <- endothelial_rna %>% 
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)

# print the Harmony PCA dimensions
ProjectDim(endothelial_rna, reduction = "harmony", dims.print = 1:5) 
# Explore heatmap of PCs
DimHeatmap(endothelial_rna, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

DimPlot(endothelial_rna, group.by = 'cell_subtype') +
  DimPlot(endothelial_rna, group.by = 'tumor_or_normal') + 
  DimPlot(endothelial_rna, group.by = 'Patient') +
  DimPlot(endothelial_rna, group.by = 'Phase') +
  DimPlot(endothelial_rna, group.by = 'seurat_clusters')

endo_markers <- FindAllMarkers(endothelial_rna)
write_csv(endo_markers, "Results/endothelial_cluster_markers_seurat.csv")

endothelial_cluster_umap <- DimPlot(endothelial_rna, group.by = 'seurat_clusters',
                        cols = DiscretePalette(n = length(
                          unique(endothelial_rna$seurat_clusters))))

endo_markers <- read_csv("Results/endothelial_cluster_markers_seurat.csv")
endothelial_markers_filtered <- endo_markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(n=5, order_by= avg_logFC)
# heatmap of the top 5 marker genes of each cluster
endothelial_cluster_heatmap <- DoHeatmap(endothelial_rna, features=endothelial_markers_filtered$gene,
                             group.colors =  DiscretePalette(n = length(
                               unique(endothelial_rna$seurat_clusters))))

endothelial_cluster_umap + endothelial_cluster_heatmap + plot_annotation(title = 'Endothelial cells')
