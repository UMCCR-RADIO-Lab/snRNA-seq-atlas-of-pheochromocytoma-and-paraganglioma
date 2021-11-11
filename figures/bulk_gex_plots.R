rm(list = ls())

library(tidyverse)
source("singlecell_colour_palettes.R")

# this version of the seurat object was created in "angio_tme_analysis.R"
# pcpg_rna <- readRDS("Data/pcpg_with_metadata.RDS")

# ----
# read in bulk RNA-seq data and metadata
# ----

# batch normalised expression
bulk_norm_expression <- read_tsv("Data/bulk_expression_batch_normalised.tsv")

# bulk metadata
bulk_metadata_all <- read_csv("Data/Table S3. Pheo-atlas metadata.csv")%>%
  mutate(Platform = "Microarray")%>%
  mutate(Platform = replace(Platform, Batch %in% c("TCGA", "Flynn_filtered"),
                            "RNA-Seq"))

bulk_metadata <- bulk_metadata_all %>%
  dplyr::select(Sample_raw, Genotype, Subtype, Malignancy) %>% 
  mutate(Subtype = recode(Subtype,!!!subtype_key2)) %>% 
  dplyr::rename("Sample" = Sample_raw)

# ----
# angio marker genes stacked violin plots
# ----

# pro-angiogenic cell type marker genes 
angio_markers <- c("DLL4", "MCAM", "VEGFA", "EPAS1", "ANGPT2", "HEY1")

angio_bulk_norm_expression <- bulk_norm_expression %>%
  filter(Gene %in% angio_markers) %>% 
  mutate(Gene = factor(Gene, levels = angio_markers))

# organise the data for plotting
angio_bulk_norm_expression_plot <- angio_bulk_norm_expression %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Count")
plot_data <-  bulk_metadata %>%
  left_join(angio_bulk_norm_expression_plot,
            by = "Sample") %>% 
  filter(Subtype != "C2C") %>% 
  mutate(Subtype = factor(Subtype, levels = subtypes_genotypes))

# plot expression by tumor subtype 
vln_angio_markers_subtypes <- ggplot(data = plot_data , aes(x = Subtype, y  = Count))+ 
  geom_violin(aes(fill = Subtype),lwd = rel(0.25)) +
  facet_wrap(~Gene,ncol = 1,
             scales = "free_y",
             strip.position="left") +
  geom_boxplot(width = rel(0.2),
               lwd = rel(0.25),
               outlier.size = rel(0.5)) +
  scale_y_continuous(position = "right",  breaks = extended_breaks(n=4)) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = subtype_genotype_cols) +
  guides(fill = "none") +
  labs(x = "", 
       y = "Normalised expression") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), 
        strip.text.y = element_text(angle = 90),
        aspect.ratio = 0.15,
        panel.spacing = unit(-0.1, "mm"),
        panel.grid=element_blank(),
        strip.background = element_blank())

vln_angio_markers_subtypes

# ggsave(vln_angio_markers_subtypes,
#        width = 150, height = 150, units = "mm",
#        filename = "Figures/angio_markers_vln_stacked.pdf")

# ----
# immune marker gene stacked violin plots
# ----

# immune cell type marker genes 
myeloid_markers <- c("CD163", "MARCO", "PLAU", "CXCL2")
lymphoid_markers <- c("CD4", "CD8A", "GNLY", "GZMB")
immune_markers <- c(myeloid_markers,   lymphoid_markers)

immune_bulk_norm_expression <- bulk_norm_expression %>%
  filter(Gene %in% immune_markers) %>% 
  mutate(Gene = factor(Gene, levels = immune_markers))

# organise the data for plotting
immune_bulk_norm_expression_plot <- immune_bulk_norm_expression %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Count")
plot_data <-  bulk_metadata %>%
  left_join(immune_bulk_norm_expression_plot,
            by = "Sample") %>% 
  filter(Subtype != "C2C") %>% 
  mutate(Subtype = factor(Subtype, levels = subtypes_genotypes))

# plot expression by tumor subtype 
vln_immune_markers_subtypes <- ggplot(data = plot_data , aes(x = Subtype, y  = Count))+ 
  geom_violin(aes(fill = Subtype),lwd = rel(0.25)) +
  geom_boxplot(width = rel(0.2),
               lwd = rel(0.25),
               outlier.size = rel(0.5)) +
  #coord_flip()+
  scale_y_continuous(position = "right", breaks = extended_breaks(n=4)) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = subtype_genotype_cols) +
  guides(fill = "none") +
  labs(x = "", 
       y = "Normalised expression") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), 
        strip.text.y = element_text(angle = 90),
        aspect.ratio = 0.15,
        panel.spacing = unit(-0.1, "mm"),
        panel.grid=element_blank(),
        strip.background = element_blank())  + 
  facet_wrap(~Gene, ncol = 1,
             scales = "free_y",
             strip.position="left") 

vln_immune_markers_subtypes

# ggsave(vln_immune_markers_subtypes,
#        width = 150, height = 150, units = "mm",
#        filename = "Figures/immune_markers_vln_stacked.pdf")

# ----
# SCLC marker gene plots
# ----

# SCLC cell type marker genes 
sclc_markers <- c("CDH19", "SOX10")

sclc_bulk_norm_expression <- bulk_norm_expression %>%
  filter(Gene %in% sclc_markers) %>% 
  mutate(Gene = factor(Gene, levels = sclc_markers))

# organise the data for plotting
sclc_bulk_norm_expression_plot <- sclc_bulk_norm_expression %>% 
  pivot_longer(cols = -Gene,
               names_to = "Sample",
               values_to = "Count")
plot_data <-  bulk_metadata %>%
  left_join(sclc_bulk_norm_expression_plot,
            by = "Sample") %>% 
  filter(Subtype != "C2C") %>% 
  mutate(Subtype = factor(Subtype, levels = subtypes_genotypes))

# plot expression by tumor subtype 

# plot expression by tumor subtype 
vln_sclc_markers_subtypes <- ggplot(data = plot_data , aes(x = Subtype, y  = Count))+ 
  geom_violin(aes(fill = Subtype),lwd = rel(0.25))+
  facet_wrap(~Gene,ncol = 6, scales = "free_x")+
  geom_boxplot(width = rel(0.2),
               lwd = rel(0.25),
               outlier.size = rel(0.5))+
  theme_bw() +
  theme(axis.text.x = element_text(size = 6), 
        aspect.ratio = 1,
        panel.spacing = unit(5, "mm"),
        text = element_text(size = 6),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6), axis.text.y = element_text(size = 6),
        panel.grid=element_blank(),
        strip.background = element_blank()) +
  coord_flip()+
  scale_fill_manual(values = subtype_genotype_cols)+
  guides(fill = "none")+
  labs(x = "Subtype", y = "Normalised expression")

vln_sclc_markers_subtypes

# ggsave(vln_SCLC_markers_subtypes,
#        width = 150, height = 100, units = "mm",
#        filename = "Figures/SCLC_markers_vln.pdf")