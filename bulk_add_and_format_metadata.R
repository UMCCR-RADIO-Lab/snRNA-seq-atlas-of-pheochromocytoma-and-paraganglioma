rm(list=ls())
library(tidyverse)
source("singlecell_colour_palettes.R")

# reformat the bulk metadata sheet so that it has the updated subtype nomenclature 
bulk_metadata <- read_csv("Data/Table S3. Pheo-atlas metadata.csv") # think i got this from the google drive

bulk_metadata <- bulk_metadata %>%
  mutate(tumor_subtype = recode(Subtype, !!!subtype_key2))
bulk_metadata$tumor_subtype <- factor(bulk_metadata$tumor_subtype, levels = c(subtypes_genotypes))

write_csv(bulk_metadata, "Data/bulk_metadata.tsv")