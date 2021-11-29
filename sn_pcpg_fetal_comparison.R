rm(list=ls())

library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)
library(ggsci)
library(scales)
source("Analysis/load_signature_genesets_jansky.R")
source("singlecell_colour_palettes.R")
subtypes_genotypes <- c("CCs (NAM)", "C1A1 (SDHx)", "C1A2 (SDHx-HN)", "C1B1 (VHL)", "C1B2 (EPAS1)", "C2A (Kinase)", "C2B1 (MAX)", "C2B2 (MAML3)", "C2C")

#####
# Analysis
#####

# ----
# Classify the neoplastic PCPG dataset 
# ----

# read in all the snRNA-seq data 
pcpg_rna <- readRDS("Data/sn.PPGLs.filtered.with.decontX.RDS")

# add the updated sample metadata 
pb_metadata <- read_tsv("Data/pseudobulk_metadata.tsv") %>% 
  mutate(Genotype_Subgroup4 = paste(Genotype, Subgroup4, sep = "_")) # make a metadata column to aggregate samples with same genotype and subgroup

# check for duplicated column names
dup.cols <- names(pb_metadata)[names(pb_metadata) %in% names(pcpg_rna[[]])]
# merge the metadata tables
metadata_new <- pcpg_rna[[]] %>%
  dplyr::select(!dup.cols[2:6]) %>% # remove the old annotations (except sample names) and replace with new ones
  left_join(pb_metadata, by = "orig.ident")
rownames(metadata_new) <- rownames(pcpg_rna[[]])
pcpg_rna@meta.data <- metadata_new

# make a subset of the PCPG RNA-seq dataset, containing neoplastic and normal adult chromaffin cells  
pcpg_rna <- subset(pcpg_rna, Cell_Type %in% c("Tumour", "Chromaffin cells", "Sustentacular cells"))
pcpg_rna[["decontX"]] <- NULL # remove this assay for extra memory
DefaultAssay(pcpg_rna) <- 'RNA'
# stash the umap coordinates in the metadata
umap_embeddings <- data.frame(Embeddings(pcpg_rna, reduction = 'umap'))
# table(rownames(umap_embeddings) == colnames(pcpg_rna))
pcpg_rna$sct_umap1 <- umap_embeddings$UMAP_1
pcpg_rna$sct_umap2 <- umap_embeddings$UMAP_2

# ----
# Re-run the seurat pipeline on subsetted data 
# ----

pcpg_rna <- NormalizeData(pcpg_rna)
pcpg_rna <- FindVariableFeatures(pcpg_rna, selection.method = "vst", nfeatures = 2000)
pcpg_rna <- ScaleData(pcpg_rna)
pcpg_rna <- RunPCA(pcpg_rna)
pcpg_rna <- RunUMAP(pcpg_rna, dims = 1:20)

# ----
# Read in the fetal data
# ----

fetal_rna <- readRDS("Data/Jansky_2021/snRNA_adrenal_medulla_jansky.RDS")
fetal_rna[["old.ident"]] <- Idents(object = fetal_rna)
fetal_rna <- FindVariableFeatures(fetal_rna)

# ----
# Classify the PCPG into the fetal cell types
# ----

# using Magnus Zethoven's method to classify the cell types based on spearman correlation (similar to scMatch) using top 3000 variable features 
jansky <- as(fetal_rna@assays$RNA@counts, "dgCMatrix")
jansky_annot <- data.frame(Idents(fetal_rna)) %>% rownames_to_column(var = "barcode") %>% dplyr::rename("ident" = "Idents.fetal_rna.")
jansky_annot <- jansky_annot[match(colnames(jansky), jansky_annot$barcode),] #%>% column_to_rownames(var = "barcode")
jansky_annot$ident <- as.character(as.factor(jansky_annot$ident))
jansky_centroids <- do.call(cbind, lapply(unique(jansky_annot$ident[!is.na(jansky_annot$ident)]) %>%
                                            (function(x){setNames(x,x)}),
                                          function(cell_type) rowMeans(jansky[, which(jansky_annot$ident == cell_type)])))
# use the 3000 most variable genes in foetal medulla
jansky_variable_features <- jansky %>% LogNormalize() %>% CreateSeuratObject() %>% FindVariableFeatures(nfeatures = 3000) %>% VariableFeatures()
f <- (jansky_annot$ident[match(colnames(jansky_centroids), jansky_annot$ident)]) %in% levels(jansky_annot$ident)
# spearman correlation, like scMatch
cluster.cells <- rownames(pcpg_rna@meta.data[pcpg_rna@meta.data$seurat_clusters %in% seq(0, 45, 1),])
c_ppgl <- cor(as.matrix(pcpg_rna@assays$SCT@data[(rownames(pcpg_rna@assays$SCT@data) %in% jansky_variable_features),(colnames(pcpg_rna@assays$SCT@data) %in% cluster.cells)]), as.matrix(jansky_centroids[rownames(pcpg_rna@assays$SCT@data[rownames(pcpg_rna@assays$SCT@data) %in% jansky_variable_features,]),]), use='pairwise.complete.obs', method="spearman")
ppgl_ident <- data.frame(barcode = rownames(c_ppgl),
                         cells = colnames(c_ppgl)[apply(c_ppgl,1,which.max)])
# write.csv(ppgl_ident, "Results/magnus_idents_var_features.csv")

ppgl_ident <- read.csv("Results/magnus_idents_var_features.csv")
# add the classifications to the seurat obj metadata
classifications <- setNames(ppgl_ident$cells,
                            ppgl_ident$barcode)
pcpg_rna <- AddMetaData(pcpg_rna,
                        col.name = "cell_type_spearman_fetal_ref_3000variable_genes",
                        metadata = classifications)
# free up memory
rm(jansky, ppgl_ident, c_ppgl, cluster.cells, f, jansky_variable_features, jansky_centroids, jansky_annot, ppgl_ident)

# ----
# Calculate gene signature scores in each of the subtypes
# ----

# calculate single-cell signature scores for each of the clusters
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(chromaffin_genes$gene),
                           name = "module_score_chromaffin")
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(latechromaffin_genes$gene),
                           name = "module_score_late_chromaffin")
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(connecting_chromaffin_genes$gene),
                           name = "module_score_connecting_chromaffin")
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(bridge_genes$gene),
                           name = "module_score_bridge")
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(neuroblast_genes$gene),
                           name = "module_score_neuroblast")
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(lateneuroblast_genes$gene),
                           name = "module_score_late_neuroblast")
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(cycling_neuroblast_genes$gene),
                           name = "module_score_cycling_neuroblast")
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(SCP_genes$gene),
                           name = "module_score_SCP")
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(late_SCP_genes$gene),
                           name = "module_score_late_SCP")
pcpg_rna <- AddModuleScore(pcpg_rna,
                           features = list(cycling_SCP_genes$gene),
                           name = "module_score_cycling_SCP")
# TODO: saveRDS(pcpg_rna, ...)

#####
# Plotting
#####

# ----
# Reformat the cell names for plotting
# ----

# fetal data cell name change
jansky_names = c("Bridge",
                 "connecting Chromaffin cells",
                 "Chromaffin cells",
                 "late Chromaffin cells",
                 "cycling Neuroblasts",
                 "Neuroblasts",
                 "late Neuroblasts",
                 "SCPs",
                 "cycling SCPs",
                 "late SCPs")

new_names <- c("Bridge",
               "Connecting progenitor cells",
               "Early chromaffin cells",
               "Late chromaffin cells",
               "Cycling neuroblasts",
               "Early neuroblasts",
               "Late neuroblasts", 
               "Early SCPs", 
               "Cycling SCPs",
               "Late SCPs") 

jansky_name_key <- setNames(new_names,jansky_names)
fetal_rna$old.ident2 <- Idents(fetal_rna)
fetal_rna$cell.type <- recode(fetal_rna$old.ident2,
                              !!!jansky_name_key)

# PCPG cell name change
md <- pcpg_rna@meta.data
md$barcode <- rownames(md)
md <- md %>% 
  mutate(Subgroup4 = case_when(
    Cell_Type == "Sustentacular cells" ~ "SCLCs", 
    Cell_Type == "Chromaffin cells" ~ "Chromaffin cells", 
    Cell_Type == "Tumour" ~ Subgroup4)) %>% 
  mutate(Sample_plot = case_when(
    Cell_Type == "Sustentacular cells" ~ "SCLCs", 
    Cell_Type == "Chromaffin cells" ~ "Chromaffin cells", 
    TRUE ~ Sample)) %>% 
  mutate(genotype_plot = case_when(
    Cell_Type == "Sustentacular cells" ~ "WT", 
    Cell_Type == "Chromaffin cells" ~ "WT", 
    TRUE ~ Genotype)) %>% 
  mutate(cell_type_spearman_fetal_ref_3000variable_genes = recode(
    cell_type_spearman_fetal_ref_3000variable_genes, !!!jansky_name_key))

pcpg_rna@meta.data <- md

pcpg_rna[["Subgroup4"]] <- pcpg_rna[["Subgroup4"]] %>%
  dplyr::mutate(Subgroup4 = recode(Subgroup4, !!!subtype_key1))

#----
# Organise colour maps  
#----

jansky.cell.cols <- setNames(object = pal_d3("category20")(20)[11:20],
                             nm = new_names)

# ----
# UMAP plots 
# ----

# ggplot theme setup
theme_set(theme_bw() + theme(panel.grid = element_blank())) # plain theme
# makes legend spacing smaller and removes legend box
t <- theme(legend.background = element_blank(),
           legend.spacing = unit(0.1, "mm"),
           legend.key.size = unit(5, "mm"), 
           plot.title = element_text(hjust = 0.5, face = "bold"))
# Makes colour legend big and square
t2 <- guides(colour = guide_legend(override.aes = list(shape = 15, size = 3)))

# enforces square ggplot panel in real units
square.ratio <- function(x){
  range2 <- function(x){sum(c(-1,1)*range(x, na.rm=T))}
  return(coord_fixed(clip="off", ratio = range2(x[,1])/range2(x[,2])))
}

sn.umap.ratio_fetal <- square.ratio(Embeddings(fetal_rna, reduction="umap"))
pt.size <- 0.15
pt.stroke <- 0.1

# umaps for the figure 
fetal_table <- cbind(fetal_rna@meta.data, Embeddings(fetal_rna, reduction='umap')) %>% 
  dplyr::rename("UMAP 1" = UMAP_1, 
                "UMAP 2" = UMAP_2)
fetal_umap <- ggplot() +
  ggtitle("Fetal adrenal medulla\n(Jansky et al.)") + 
  geom_point(data = fetal_table, mapping = aes(x = `UMAP 1`, y = `UMAP 2`, col = cell.type), size = pt.size, stroke = pt.stroke) +
  t + 
  t2 +
  scale_color_manual(values = jansky.cell.cols) +
  sn.umap.ratio_fetal +
  labs(colour="Cell type") + 
  theme(legend.text.align = 0,
        legend.position = "bottom")

sn.umap.ratio_pcpg <- square.ratio(Embeddings(pcpg_rna, reduction="umap"))
pcpg_table <- cbind(pcpg_rna@meta.data, Embeddings(pcpg_rna, reduction='umap')) %>% 
  dplyr::rename("UMAP 1" = UMAP_1, 
                "UMAP 2" = UMAP_2)

# arrange the table so that the rarer classifications are plotted on top
pcpg_table <- pcpg_table %>%
  group_by(cell_type_spearman_fetal_ref_3000variable_genes) %>%
  mutate(cell_type_count=n()) %>%
  ungroup() %>%
  arrange(-cell_type_count)

pcpg_umap <- ggplot() +
  ggtitle("PCPG") + 
  geom_point(data = pcpg_table, mapping = aes(x = `UMAP 1`, y = `UMAP 2`, col = cell_type_spearman_fetal_ref_3000variable_genes), size = pt.size, stroke = pt.stroke) +
  t + 
  t2 +
  scale_color_manual(values = jansky.cell.cols) +
  sn.umap.ratio_pcpg +
  labs(colour="Predicted cell type") + 
  theme(legend.position = "none", legend.text.align = 0)

combined_umaps <- fetal_umap + pcpg_umap +
  plot_layout(guides = "collect") & theme(plot.title = element_text(hjust = 0.5),
                                          text = element_text(size=6),
                                          legend.title = element_text(size=6),
                                          legend.text = element_text(size =5),
                                          axis.text = element_text(size=5),
                                          legend.key.size = unit(6, "points"))
combined_umaps

# ----
# Classifications stacked bar plot
# ----

# proportions of samples 
sample.data.vargenes3000 <- pcpg_rna@meta.data %>%
  dplyr::select(cell_type_spearman_fetal_ref_3000variable_genes,
                Sample_plot,
                genotype_plot,
                Subgroup4) %>%
  group_by(Sample_plot, Subgroup4, genotype_plot, cell_type_spearman_fetal_ref_3000variable_genes) %>%
  summarise(vargenes3000_count = n()) %>% # cell count in groups
  mutate(total = sum(vargenes3000_count)) %>%
  mutate(vargenes3000_prop = vargenes3000_count / total)

# faceted sample plots ----

sample_vargenes3000_pct <- ggplot(sample.data.vargenes3000,
                                  aes(x = Sample_plot,
                                      y = vargenes3000_prop,
                                      fill = cell_type_spearman_fetal_ref_3000variable_genes)) +
  geom_col(position = 'fill') +
  theme_classic() +
  facet_grid(~Subgroup4 + genotype_plot, scales = "free_x", space = "free_x")+
  labs(y = "Fraction") +
  scale_fill_manual(values = jansky.cell.cols) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        strip.background = element_blank(),
        panel.grid = element_blank())

sample_vargenes3000_pct +
  labs(colour = "") +
  plot_layout(guides = 'collect') 

# ggsave(sample_vargenes3000_pct, 
#   width = 300, height = 100, units = "mm",
#   filename = "Figures/stacked_bar_fetal_label_transfer_magnusmethod.pdf",
#   useDingbats = FALSE )

# ----
# module scores stacked violin plot 
# ----

# the features I want to plot are: 
sn_features <- c("module_score_bridge1",
                 "module_score_connecting_chromaffin1",
                 "module_score_chromaffin1",
                 "module_score_late_chromaffin1",
                 "module_score_neuroblast1",
                 "module_score_late_neuroblast1",
                 "module_score_SCP1",
                 "module_score_late_SCP1",
                 "module_score_cycling_SCP1")

cell_types <- c("Bridge",
                "Connecting progenitor cells",
                "Early chromaffin cells",
                "Late chromaffin cells",
                "Early neuroblasts",
                "Late neuroblasts",
                "Early SCPs",
                "Late SCPs",
                "Cycling SCPs")
# key for converting the module names back to the relevant cell type
cell_types_key <- setNames(cell_types, sn_features)

# get the plot data from seurat object metadata 
sn_data <- pcpg_rna[[c("Sample",
                       "orig.ident",
                       "Cell_Type",
                       "Subgroup4",
                       "Genotype",
                       sn_features)]]
# stash the barcodes
sn_data$barcode <- rownames(sn_data)

# make a metadata column for plotting 
sn_data <- sn_data %>%
  tibble() %>% 
  mutate(Subgroup4 = case_when(
    Cell_Type == "Sustentacular cells" ~ "SCLCs", 
    Cell_Type == "Chromaffin cells" ~ "Chromaffin cells", 
    Cell_Type == "Tumour" ~ recode(Subgroup4, !!!subtype_key2))) %>%
  mutate(Subgroup4 = factor(Subgroup4,
                            levels = rev(c("Chromaffin cells",
                                           "SCLCs",
                                           subtype_key1)))) %>% 
  # compress module scores into a single column, with the module_cell_type column designating the cell type
  pivot_longer(names_to = "module_cell_type",
               values_to = "module_score",
               cols = starts_with("module_score_")) %>% 
  mutate(module_cell_type = recode(module_cell_type, !!!cell_types_key)) %>% 
  mutate(module_cell_type = factor(module_cell_type, levels = jansky_name_key))

sample.order <- sn_data %>%
  arrange(Subgroup4, Genotype, Sample) %>%
  dplyr::select(Sample) %>%
  distinct() %>%
  pull(Sample)

# Draw plots
col_palette <- c(subtype_genotype_cols, "Chromaffin cells"="light pink", "SCLCs"="light pink") 
# plot expression by tumor subtype 
vln_module_scores <- ggplot(data = sn_data , aes(x = Subgroup4, y = module_score))+ 
  geom_violin(aes(fill = Subgroup4),lwd = rel(0.25)) +
  facet_wrap(~module_cell_type,ncol = 9,
             scales = "free_x",
             strip.position="top") +
  coord_flip()+
  scale_y_continuous(position = "left",limits = c(-0.3,1.3), breaks = extended_breaks(n=3)) +
  # scale_x_discrete(position = "bottom") +
  scale_fill_manual(values = col_palette) +
  guides(fill = "none") +
  labs(x = "", 
       y = "Normalised expression") +
  theme_bw() +
  theme(
    #axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), 
    strip.text = element_text(angle = 90, hjust = 0, vjust =0.5),
    aspect.ratio = 5,
    panel.spacing = unit(-0.1, "mm"),
    panel.grid=element_blank(),
    strip.background = element_blank())

vln_module_scores
