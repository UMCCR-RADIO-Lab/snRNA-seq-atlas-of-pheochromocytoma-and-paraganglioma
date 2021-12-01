#setwd("/homevol/admin/snPPGL/")
#.libPaths(c("Rlibs"))

source("/researchers/magnus.zethoven/Seurat/seurat3.6.R")
library(dplyr)
library(ggplot2)
library(Matrix)

z_score = function(x, center=TRUE, scale=TRUE, robust=FALSE, log=FALSE){
  if(log){
    x = log(x)
  }
  if(robust){
    u = median(x)
    a = mad(x)
    if(a == 0){
      a = 1.253314*mean(abs(x - median(x)))
    }
    return((x - ifelse(center, u, 0))/(ifelse(scale, a, 0)))
  } else {
    return(scale(x, center, scale))
  }
  
}


library(doParallel)
c25 = c("dodgerblue2","#E31A1C", # red
        "green4",
        "#6A3D9A", # purple
        "#FF7F00", # orange
        "gold1",
        "skyblue2","#FB9A99", # lt pink
        "palegreen2",
        "#CAB2D6", # lt purple
        "#FDBF6F", # lt orange
        "gray70", "khaki2",
        "maroon","orchid1","deeppink1","blue1","steelblue4",
        "darkturquoise","green1","yellow4","yellow3",
        "darkorange4","brown",
        "black")

c30 = c(c25, c("cyan", "lightsteelblue", "goldenrod4", "navyblue", "green4"))
library(future)
Sys.setenv(R_LIBS = paste(.libPaths()[1], Sys.getenv("R_LIBS"), sep=.Platform$path.sep))
plan("sequential")
options(future.globals.maxSize = 24000 * 1024^2)


# load raw data
load("data/sn.PPGLs.Rdata")
samples = unique(sn.PPGLs$orig.ident)
annot = read.delim("data/PPGL_Single_cell_metadata.tsv", sep='\t', header=T, as.is=T)
meta.data = read.delim('data/meta.data2.tsv', sep='\t', header=T, as.is=T,row.names=1) 
meta.data$Barcode = rownames(meta.data)
sn.PPGLs = sn.PPGLs %>% PercentageFeatureSet(pattern="^MT-", col.name = "percent.mt")

sn.PPGLs@meta.data[,c("Scrublet", "scMatch", "Cell_Type_Old")] = meta.data[,c("Scrublet", "scMatch", "Assumed_Cell_Type")]
x = sn.PPGLs@meta.data
x = x %>% group_by(orig.ident) %>% mutate(nCount_RNA_MAD = z_score(nCount_RNA, log=TRUE, robust=TRUE),
                                          nFeature_RNA_MAD = z_score(nFeature_RNA, log=TRUE, robust=TRUE),
                                          percent.mt_MAD = z_score(percent.mt, robust=TRUE),
                                          Scrublet_MAD = z_score(Scrublet, robust=TRUE))
sn.PPGLs@meta.data[,c("nCount_RNA_MAD", "nFeature_RNA_MAD", "percent.mt_MAD", "Scrublet_MAD")] = x[,c("nCount_RNA_MAD", "nFeature_RNA_MAD", "percent.mt_MAD", "Scrublet_MAD")]


# Filter low quality cells and doublets
samples.to.remove = "A3"
# Oops, forgot neutrophils >_<
mad.thresh = ifelse(grepl("B cell|T cell|mast cell|natural killer", sn.PPGLs$scMatch), -4, -2.5)
sn.PPGLs = sn.PPGLs %>% subset(cells = colnames(sn.PPGLs)[which(sn.PPGLs$percent.mt_MAD <= 5 & sn.PPGLs$Scrublet_MAD <= 2 & sn.PPGLs$nCount_RNA_MAD >= mad.thresh & sn.PPGLs$nFeature_RNA_MAD >= mad.thresh & !(sn.PPGLs$orig.ident %in% samples.to.remove))])

sn.PPGLs = sn.PPGLs %>% SCTransform(verbose = TRUE, variable.features.n = NULL, conserve.memory = T, latent_var_nonreg = "percent.mt")
sn.PPGLs = sn.PPGLs %>% CellCycleScoring(s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
sn.PPGLs = sn.PPGLs %>% SCTransform(verbose = TRUE, variable.features.n = NULL, conserve.memory = T, latent_var_nonreg = c("percent.mt", "G2M.Score", "S.Score"))
sn.PPGLs = sn.PPGLs %>% RunPCA(npcs=80) 

metric="cosine"
dims = 1:20
k=20
resolution=0.8

sn.PPGLs = sn.PPGLs %>%
  FindNeighbors(reduction="pca", dims=dims, nn.method="annoy", annoy.metric = metric, k.param = k) %>%
  RunUMAP(dims=dims, metric=metric, n.neighbors = k) %>%
  FindClusters(resolution=resolution)

saveRDS(sn.PPGLs, file = "data/sn.PPGLs.sct.latent_var_nonreg.rds")

sn.PPGLs$Cell_Type_Old = gsub(" \\(.*$", "", gsub("Adipocyte", "Fibroblast", sn.PPGLs$Cell_Type_Old))
t = table(sn.PPGLs$Cell_Type_Old, Idents(sn.PPGLs))
sn.PPGLs$Cell_Type = apply(t, 2, function(x){rownames(t)[which(x == max(x))[1]]})[match(Idents(sn.PPGLs), colnames(t))]
rm(t)
sn.PPGLs$Cell_Type[sn.PPGLs$Cell_Type == "Tumour" & grepl("^NAM", sn.PPGLs$orig.ident)] = "Chromaffin cells"

saveRDS(sn.PPGLs, file = "data/sn.PPGLs.sct.latent_var_nonreg.rds")
