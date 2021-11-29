# NATMI analysis + figures
# Natmi Plots - Sustentacular cells and subtype-specific receptors volcano plot

library(tidyverse)
library(circlize)
library(ggrepel)
library(patchwork)
library(RColorBrewer)
library(scales)
source("Figures/dotplot_functions.R")
source("singlecell_colour_palettes.R")

# Ggplot2 blank theme
blank_theme <- theme_bw()+ #base_size = 25
  theme(panel.grid=element_blank(),
        strip.background = element_blank())

# ----
# Read in data 
# ----

# Read in the NATMI results:
All_natmi <- read_csv("Data/NATMI/Table S11 Per-sample NATMI results.csv")
# Read in the pseudobulk metadata
sn_md <- read_tsv("Data/pseudobulk_metadata.tsv")

# Get the genes upregulated in SCLCs compared to all other cells
# Pval < 0.05 and logFC > 3

all_dge <- read_csv("Data/Pseudobulk all DGE.csv")

sust_specific_genes <- all_dge %>%
  filter(Contrast == "Sustentacular.cells_normal" &
           adj.P.Val < 0.05 &
           logFC > 3) %>%
  arrange(-logFC) 


# ----
# preprocess and summarise NATMI results:
# ----

# filter natmi results: 
#  ligand and receptor expressed in at least 10 cells 
#  Ligand and receptor detection rate over 10% 
clean_natmi <- All_natmi %>%
  filter(`Ligand-expressing cells` > 10 &
           `Receptor-expressing cells` > 10 &
           `Ligand detection rate` > 0.1 & `Receptor detection rate` > 0.1 &
           `Sending cluster`!=`Target cluster` &
           `Ligand symbol` != `Receptor symbol`) %>%
  mutate(T_N = ifelse(grepl("NAM", Sample), "Normal", "Tumour"))

# write_tsv(clean_natmi, "Results/natmi_table_clean.tsv")

# ----
# Plot the differentially expressed receptors
# ----

receptors <- read_tsv("Data/NATMI/Receptor_database.txt") %>%
  # Keep only receptors
  filter(grepl("Receptor", Type)) %>%
  filter(Type != "Receptor/Ligand")

tumour_specific_genes <- all_dge %>%
  mutate(Contrast = recode(Contrast, !!!subtype_key1)) %>% 
  filter(Contrast %in% c("CCs (NAM)", "C1A1 (SDHx)", "C1A2 (SDHx-HN)", "C1B1 (VHL)", "C1B2 (EPAS1)", "C2A (Kinase)", "C2B1 (MAX)", "C2B2 (MAML3)")) %>% 
  filter(adj.P.Val < 0.05, logFC > 0) %>%
  arrange(-logFC) %>% 
  mutate(is_receptor = if_else(Gene %in% receptors$`Hgnc Symbol`, TRUE, FALSE)) %>% 
  mutate(significant = if_else(adj.P.Val < 0.05 & logFC > 3, TRUE, FALSE)) 

# make a volcano plot 
# receptors are coloured by genotype
# non significant genes

tumour_receptors_volcano <- ggplot(tumour_specific_genes,
                                   aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(data = filter(tumour_specific_genes, is_receptor == FALSE),  col = "lightgrey") +
  geom_point(data = filter(tumour_specific_genes, is_receptor == TRUE & significant == TRUE),  aes(col = Contrast)) +
  geom_text_repel(data = tumour_specific_genes %>%
                    filter(is_receptor == TRUE & significant == TRUE) %>%
                    group_by(Contrast) %>% 
                    slice_max(order_by = logFC, n = 5) %>% 
                    ungroup(),
                  aes(colour = Contrast, label = Gene),
                  size = 6,
                  max.overlaps = Inf, box.padding = 0.55, show.legend = FALSE) + 
  geom_vline(xintercept = 3, linetype = 'dashed', col = 'red') + 
  theme_bw() +
  labs(x = "Log2FC") +
  theme(panel.grid=element_blank(),
        strip.background = element_blank(),
        text=element_text(size=16),
        axis.text = element_text(size=16),
        aspect.ratio=1) + 
  scale_colour_manual(values = subtype_genotype_cols[2:8])

tumour_receptors_volcano
# 
# ggsave(tumour_receptors_volcano,
#        filename = "Figures/tumour_receptors_volcano.pdf",
#        height = unit(7.5, "cm"),
#        width = unit(10, "cm"))

# ----
# Investigate sustentacular cell signalling
# ----

# get natmi interactions involving sustentacular cells 
sust_interactions_natmi <- clean_natmi %>%
  filter(`Sending cluster` == "Sustentacular cells" | `Target cluster` ==  "Sustentacular cells")

# get the natmi results for sust-specific genes
sust_tumor_interactions_natmi <- sust_interactions_natmi %>% 
  # keep only ligands / receptors that are over-expressed in sustentacular cells 
  filter((`Sending cluster` == "Sustentacular cells" & `Target cluster` == "Tumour" & `Ligand symbol` %in% sust_specific_genes$Gene ) | # & `Receptor symbol` %in% tumour_specific_genes$Gene
           (`Target cluster` ==  "Sustentacular cells" & `Sending cluster` == "Tumour" & `Receptor symbol` %in% sust_specific_genes$Gene)) %>% # & `Ligand symbol` %in% tumour_specific_genes$Gene
  dplyr::rename(orig.ident = Sample) %>%
  left_join(sn_md, by = "orig.ident") %>% # join to sample metadata
  # calculate mean expression weight for each interaction
  group_by(`Sending cluster`, `Target cluster`,`Ligand symbol`,`Receptor symbol`) %>%
  summarise(weight = mean(`Edge total expression weight`)) %>%
  ungroup() %>% 
  arrange(`Sending cluster`, -weight) %>% 
  mutate(`Ligand symbol` = factor(`Ligand symbol`, levels = unique(`Ligand symbol`))) %>%
  mutate(`Receptor symbol` = factor(`Receptor symbol`, levels = unique(`Receptor symbol`)))

natmi_plot <- ggplot(sust_tumor_interactions_natmi, aes(x=`Ligand symbol`, y = `Receptor symbol`, fill = log2(weight)))+
  geom_tile()+
  coord_equal() +
  theme_bw()+
  # Change the angle of the x axis labels
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(),
        #panel.grid=element_blank(),
        legend.position = "bottom",
        strip.background = element_blank())+
  scale_fill_gradient2(low = "white",mid = "orange", high = "red",midpoint = 33.5) +
  labs(fill = "Log2 edge total\nexpression weight", x = "Ligand symbol", y = "Receptor symbol") +
  guides(fill = guide_colourbar(barheight = 0.5, barwidth = 5, direction = "horizontal"))

natmi_plot

# ggsave(natmi_plot,
#        filename = "Figures/NATMI_SCLC_Tumour_heatmap.pdf",
#        height = unit(10, "cm"),
#        width = unit(15, "cm"))
