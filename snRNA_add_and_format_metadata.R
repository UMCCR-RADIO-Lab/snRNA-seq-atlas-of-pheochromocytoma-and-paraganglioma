rm(list = ls())

library(Seurat)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(scales)
library(circlize)
library(UpSetR)

source("Figures/dotplot_functions.R")
source("singlecell_colour_palettes.R")

# ----
# Read in the data and metadata and reformat the metadata 
# ----

# read in all the snRNA-seq data 
pcpg_rna <- readRDS("Data/sn.PPGLs.filtered.with.decontX.RDS")

# add the updated sample metadata 
pb_metadata <- read_tsv("Data/pseudobulk_metadata.tsv") %>% 
  mutate(Genotype_Subgroup4 = paste(Genotype, Subgroup4, sep = "_")) # make a metadata column to aggregate samples with same genotype and subgroup

# duplicated columns 
dup.cols <- names(pb_metadata)[names(pb_metadata) %in% names(pcpg_rna[[]])]
# merge the metadata tables
metadata_new <- pcpg_rna[[]] %>%
  dplyr::select(!dup.cols[2:6]) %>% # remove the old annotations (except sample names) and replace with new ones
  left_join(pb_metadata, by = "orig.ident") 

rownames(metadata_new) <- rownames(pcpg_rna[[]])

# recode cell type labels 
metadata_new <-  metadata_new %>%
  mutate(Cell_Type = recode(Cell_Type, 
                            "Macrophages/Monocytes" = "Myeloid cells",
                            "T cells" = "T/NK cells", 
                            "Sustentacular cells" = "SCLCs")) 

# add in the cell subtype scMatch labels 
cell_subtype_scmatch <- read_tsv("/data/gpfs/projects/punim0010/projects/Pattison_projects/PPGL_sc_RNA_seq_compendium/Metadata/pseudobulk_cell_types.tsv")
cell_subtype_scmatch$Cell_Type <- gsub("pos", "+", cell_subtype_scmatch$Cell_Type)
cell_subtype_scmatch$Cell_Type <- gsub("neg", "-", cell_subtype_scmatch$Cell_Type)

table(cell_subtype_scmatch$Barcode == rownames(metadata_new)) # the barcodes are all the same 
metadata_new$cell_subtype <- cell_subtype_scmatch$Cell_Type

# recode the cluster labels 
subtypes_old = c("Normal", "C1Ai", "C1Aii", "C1Bi", "C1Bii", "C2A", "C2Bi", "C2Bii", "C2C")
subtypes <- setNames(subtypes_genotypes, subtypes_old)
metadata_new <- metadata_new %>% 
  mutate(tumor_subtype = recode(Subgroup4, !!!subtypes)) %>%
  group_by(tumor_subtype, cell_subtype) %>%
  mutate(number_cell_subtype_in_tumor_subtype = n()) %>% 
  ungroup() %>% 
  data.frame()

rownames(metadata_new) <- rownames(pcpg_rna[[]])

pcpg_rna@meta.data <- metadata_new
DimPlot(pcpg_rna, group.by = "cell_subtype")

pcpg_rna[["decontX"]] <- NULL # remove this assay for extra memory
# stash the umap coordinates in the metadata
umap_embeddings <- data.frame(Embeddings(pcpg_rna, reduction = 'umap'))
pcpg_rna$sct_umap1 <- umap_embeddings$UMAP_1
pcpg_rna$sct_umap2 <- umap_embeddings$UMAP_2

pcpg_rna$tumor_or_normal <- if_else(pcpg_rna$Sample %in% c("E240", "E243"), "Normal", "Tumor")

# read in scores for all the barcodes - calculated in the fetal_comparison script
module_scores <- readRDS("Results/all_module_scores.csv") %>% 
  column_to_rownames(var = "barcode") %>% 
  data.frame()

pcpg_rna <- AddMetaData(pcpg_rna,
                        metadata = module_scores) 

# read in the table with batch metadata
snRNA_qc <- read_csv("Data/Table S2. snRNA-Seq QC data.csv")

md <- pcpg_rna@meta.data
md$barcode <- rownames(md)

# join to existing metadata and add to the seuratobj
md <- md %>%
  left_join(snRNA_qc) %>% 
  select(-Batch, -Subgroup2, -Subgroup, -Subgroup3, -Subgroup4, -Genotype2) %>% # remove old metadata columns
  mutate(pseudohypoxic_or_kinase_signalling_samplelevel = substr(tumor_subtype, start = 0, stop = 2)) %>%  # add a column for C1 vs C2 cluster of the genotype 
  mutate(pseudohypoxic_or_kinase_signalling = case_when(
    pseudohypoxic_or_kinase_signalling_samplelevel == "CC" ~ "Normal",
    Cell_Type != "Tumour" ~ "Normal",
    TRUE ~ pseudohypoxic_or_kinase_signalling_samplelevel)) %>% 
  mutate(Patient = str_remove(Sample, pattern = "-.+")) %>% 
  mutate(cell_subtype = recode(cell_subtype,
                                "Sustentacular cells" = "SCLCs")) %>% 
  mutate(Sample = recode(Sample,"E326" = "E236" )) # this sample name had a typo in it 

rownames(md) <- md$barcode
pcpg_rna@meta.data <- md

# Add the myeloid dim reduction that magnus did 
myeloid_umap <- readRDS("Data/myeloid_cells_umap.rds")

pcpg_rna[["Myeloid.cells"]] <- CreateDimReducObject(embeddings = myeloid_umap,
                                                    key = "Myeloid_",
                                                    assay = DefaultAssay(pcpg_rna))

lymphoid_umap <- readRDS("Data/Lymphoid_cells_umap.rds")

pcpg_rna[["Lymphoid.cells"]] <- CreateDimReducObject(embeddings = lymphoid_umap,
                                                    key = "Lymphoid_",
                                                    assay = DefaultAssay(pcpg_rna))

saveRDS(pcpg_rna, "Data/pcpg_with_metadata_and_qc.RDS")
