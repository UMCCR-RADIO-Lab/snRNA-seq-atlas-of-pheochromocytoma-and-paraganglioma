rm(list=ls())
# inferCNV plots
library(infercnv)
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(fastcluster)
library(RColorBrewer)
library(magick)

source("Figures/dotplot_functions.R")
ht_opt$fast_hclust = TRUE # use fastcluster for heatmap clustering

# constants -------------------------------------------------------------------------------------------

infercnv_rds_folder <- "Data/infercnv/infercnv_objects"

# colors
col_fun1 = colorRamp2(
  breaks = c(-.1,0,.1),
  colors=c("blue", "white","red"))

col_fun_snp6 <- colorRamp2(
  breaks = seq(-0.3, 0.3, length.out = 11),
  colors = rev(brewer.pal(11, "RdBu")))


# functions -------------------------------------------------------------------------------------------

read_infercnv_rds <- function(sample_name,
                              infercnv_dir){
  #' read in the infercnv object
  infercnv_file <- file.path(infercnv_dir, sample_name, "run.final.infercnv_obj")
  infercnv_obj <- readRDS(infercnv_file) # this has all normal cells from all samples (normal reference) plus the SCLCs and tumor cells from the test sample  
  return(infercnv_obj)
}


plot_infercnv_heatmap <- function(
  infercnv_obj,
  metadata,
  sample_name,
  height = 5,
  width = 10,
  outdir,
  save = FALSE,
  segments_file){
  #' make a complexheatmap from an infercnv object
  #' partly modified from the infercnv plot_cnv function at https://github.com/broadinstitute/infercnv/blob/master/R/inferCNV_heatmap.R
  #' @param infercnv: an infercnv object
  #' @param metadata: cell-level metadata with cell barcodes as row names 
  #' @param sample_name: name of the tumor sample, 
  #' @param outdir: directory in which to save the plot
  #' @param save: TRUE to save heatmap as a PNG
  
  # get the pseudo expression matrix from the infercnv object
  plot_data <- infercnv_obj@expr.data
  
  # automatically generate the plotting thresholds ---------------------------------------------------------------------
  # this data transformation is the same as performed by the inferCNV package plotting function
  x.center <- mean(infercnv_obj@expr.data)
  # examine distribution of data that's off-center, since much of the center could
  # correspond to a mass of data that has been wiped out during noise reduction
  quantiles = quantile(plot_data[plot_data != x.center], c(0.01, 0.99))
  # determine max distance from the center.
  delta = max( abs( c(x.center - quantiles[1],  quantiles[2] - x.center) ) )
  low_threshold = x.center - delta
  high_threshold = x.center + delta
  x.range = c(low_threshold, high_threshold)
  
  print(paste0("x.center: ", x.center))
  print(paste("auto thresholding at:", low_threshold, high_threshold, sep = " "))
  
  # replace extreme values with the threshold values
  plot_data[plot_data < low_threshold] <- low_threshold
  plot_data[plot_data > high_threshold] <- high_threshold
  
  # get the barcodes for reference cells and observation cells 
  all_barcodes <- colnames(plot_data)
  sample_barcodes <- metadata[metadata$Sample == sample_name,"barcode"]
  
  # reference cells heatmap ---------------------------------------------------------------------
  
  # get the reference cell indices 
  ref_idx <- unlist(infercnv_obj@reference_grouped_cell_indices) 
  ref_idx = ref_idx[order(ref_idx)]
  ref_barcodes <- all_barcodes[ref_idx] # reference cell barcodes
  
  # subset reference cells from the plot matrix
  ref_data <- plot_data[, ref_barcodes, drop=FALSE]
  # subset the seurat metadata to get the cell types for reference cells 
  ref_cell_types <- metadata[ref_barcodes, "Cell_Type", drop=FALSE] # for cell type anno
  
  ref_cell_type_anno <-  HeatmapAnnotation(
    df = ref_cell_types,
    col = list("Cell_Type" = cell.cols),
    which = "row",
    show_annotation_name = FALSE)
  
  ref_data <- t(ref_data - 1) 
  
  # make the reference heatmap ---------------------------------------------------------------------
  ref_hm <- Heatmap(
    ref_data,
    width = unit(width, "cm"),
    height = unit(height, "cm"),
    col = col_fun1,
    right_annotation = ref_cell_type_anno,
    name="InferCNV reference",
    cluster_columns = FALSE,
    cluster_rows=FALSE,
    show_row_names = FALSE,
    show_column_names=FALSE,
    use_raster=TRUE,
    raster_by_magick=FALSE,
    raster_quality = 2,
    row_split = ref_cell_types$Cell_Type,
    row_title_rot = 0,
    show_heatmap_legend=FALSE
  )
  
  # observation cells heatmap ---------------------------------------------------------------------
  
  # subset observation cells from the plot matrix
  obs_barcodes <- all_barcodes[-ref_idx] # observation cell barcodes
  
  # print(paste0("obs_barcodes: ", length(obs_barcodes), "\nsample_obs_barcodes: ",  length(sample_obs_barcodes)))
  obs_data <- plot_data[,obs_barcodes, drop=FALSE]
  
  # obs_data <- scale(obs_data)
  obs_data <- t(obs_data - 1) 
  
  # observation cell type annotation ---------------------------------------------------------------------
  # subset the seurat metadata to get the cell types for observation cells 
  obs_cell_types <- metadata[obs_barcodes, "Cell_Type", drop=FALSE] # for cell type anno
  
  # subset the seurat metadata to get the cell types for reference cells 
  obs_cell_types <- metadata[obs_barcodes, "Cell_Type", drop=FALSE]
  obs_cell_type_anno <-  HeatmapAnnotation(
    df = obs_cell_types,
    col = list("Cell_Type" = cell.cols),
    which = "row",
    show_annotation_name = FALSE)
  
  # observation clustering tree ---------------------------------------------------------------------------
  
  print(paste0("clustering with euclidean distance and ward.D2 method"))
  d <- dist(obs_data, method = "euclidean")
  obs_hcl <- fastcluster::hclust(d, method = "ward.D2")
  print(obs_hcl)
  
  # observation clustering annotation bar ---------------------------------------------------------------------
  k <- if_else(sum(obs_cell_types == "SCLCs") == 0, 2, 3)
  subclusters <- data.frame(subclusters=as.character(cutree(obs_hcl, k=k)))
  subcluster_colours <- setNames(object = brewer.pal(name = "Dark2", n = 3),
                                 nm = as.character(1:3))
  print(subcluster_colours)
  
  subcluster_anno <- HeatmapAnnotation(
    df = subclusters,
    name='subclusters',
    which="row",
    col = list('subclusters' = subcluster_colours))
  
  # make the observation heatmap ---------------------------------------------------------------------
  
  obs_hm <- Heatmap(
    obs_data,
    width = unit(width, "cm"),
    height = unit(height, "cm"),
    right_annotation = obs_cell_type_anno,
    left_annotation = subcluster_anno,
    name = paste0(sample_name, " inferCNV"),
    col = col_fun1,
    cluster_columns = FALSE,
    cluster_rows = obs_hcl,
    cluster_row_slices = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    use_raster = TRUE,
    raster_by_magick = FALSE,
    raster_quality = 2,
    raster_device="CairoPNG",
    row_title_rot = 0
  )
  
  # make an annotation for the chromosomes ---------------------------------------------------------------------

  contigs <- infercnv_obj@gene_order[["chr"]]  # Contigs
  unique_contigs <- unique(contigs)
  n_contig <- length(unique_contigs)
  ct.colors <- setNames(object = rep(brewer.pal(n=3, "Accent")[2:3],
                                     n_contig/2),
                        nm = unique_contigs)
  
  # use this if no snp6 array data is available
  # chrom_anno <- HeatmapAnnotation(
  #   df = data.frame(
  #     Chromosome = contigs
  #   ),
  #   col = list(
  #     Chromosome = ct.colors
  #   ),
  #   show_legend = FALSE
  # )
  # ht_list <- ref_hm %v% obs_hm %v% chrom_anno
  
  # if the snp6 data is specified then plot it
  
  # make an annotation for the SNP6 array data ----------------------------------------------------
  segments_table <- read.delim(segments_file,
                               header=T,
                               as.is=T,
                               sep='\t')
  
  genes <- rownames(infercnv_obj@expr.data)
  snp6_array <- sapply(genes, function(i){
    a=infercnv_obj@gene_order[i,,drop=F]
    b=segments_table[which(segments_table$Chromosome==a[,1]),,drop=F]
    return(tail(b$Value[b$Start < a$stop], 1))
  })
  
  # combined annotation of snp6 and chromosomes ----------------------------------------------------
  
  bottom_anno_df <- data.frame(
    SNP6 = snp6_array,
    Chromosome = contigs
    )
  
  bottom_anno <- HeatmapAnnotation(
    df = bottom_anno_df,
    col = list(
      SNP6 = col_fun_snp6,
      Chromosome = ct.colors
    ),
    show_legend = c(TRUE,FALSE))
  
  # combine heatmap components
  ht_list <- ref_hm %v% obs_hm %v% bottom_anno
  
  # draw heatmaps ---------------------------------------------------------------------
  
  if (save == "PNG"){
    outfile <-  file.path(
      outdir,
      paste0("infercnv_heatmap_", sample_name, ".png"))
    print(paste0("Saving at: ", outfile))
    png(outfile, width=20, height=20, units="cm", res = 1000)
    draw(ht_list,  merge_legends=TRUE)
    dev.off()
  }else if (save){
    outfile <- file.path(
      outdir,
      paste0("infercnv_heatmap_", sample_name, ".pdf"))
    print(paste0("Saving at: ", outfile))
    pdf(outfile,
        width=20,
        height=15)
    draw(ht_list, merge_legends=TRUE)
    dev.off()
  }
}

# read in files -------------------------------------------------------------------------------------------

# get singlecell metadata
pcpg_rna <- readRDS("Data/pcpg_with_metadata.RDS") 
md <- pcpg_rna@meta.data
md$barcode <- rownames(md)

md <- md %>%
  select(-Batch, -Subgroup2, -Subgroup, -Subgroup3, -Subgroup4, -Genotype2) %>% # remove old metadata columns
  mutate(pseudohypoxic_or_kinase_signalling = substr(tumor_subtype, start = 0, stop = 2)) %>%  # add a column for C1 vs C2 cluster of the genotype 
  mutate(Patient = str_remove(Sample, pattern = "-.+")) 
rownames(md) <- md$barcode
rm(pcpg_rna)

# these were generated by running Magnus' infercnv script
PGL1_infercnv <- read_infercnv_rds("PGL1",
                                   infercnv_rds_folder)
PGL3_infercnv <- read_infercnv_rds("PGL3",
                                   infercnv_rds_folder)
E215_infercnv <- read_infercnv_rds("A11",
                                   infercnv_rds_folder)
E024_infercnv <- read_infercnv_rds("VHP35T",
                                   infercnv_rds_folder)
E042_infercnv <- read_infercnv_rds("VCBP14T",
                                   infercnv_rds_folder)

# samples with SNP6 
plot_infercnv_heatmap(infercnv_obj=PGL1_infercnv,
                      metadata=md,
                      sample_name="P018-PGL1",
                      outdir = "Figures",
                      save=TRUE, 
                      segments_file = "Data/SNP6_segments/10_119_P018_P1-1_(CytoScanHD_Array)/segments.txt")

plot_infercnv_heatmap(infercnv_obj=E024_infercnv,
                      metadata=md,
                      sample_name="E024",
                      outdir="Figures",
                      save=TRUE,
                      segments_file = "Data/SNP6_segments/AF_VPH35T_(CytoScanHD_Array)/segments.txt")

plot_infercnv_heatmap(infercnv_obj=E042_infercnv,
                      metadata=md,
                      sample_name="E042",
                      outdir="Figures",
                      save=TRUE,
                      segments_file = "Data/SNP6_segments/AF_VCBP14T_(CytoScanHD_Array)/segments.txt")

plot_infercnv_heatmap(infercnv_obj=PGL3_infercnv,
                      metadata=md,
                      sample_name="P018-PGL3",
                      outdir = "Figures",
                      save=TRUE,
                      segments_file = "Data/SNP6_segments/GSE94378_CytoScanHD_segments/GSM2474324_C1_PGL3_CytoScanHD_Array__segments.txt")
