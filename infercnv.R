source("/researchers/magnus.zethoven/Seurat/seurat3.6.R")
#load("data/sn.PPGLs.Rdata")
#source("funcs.R")
plan("multiprocess", workers = 1)
options(future.globals.maxSize = 20 * 1000 * 1024^2)


sample = commandArgs(trailingOnly = T)[1]
print(sample)

normal.cell.types = c("Adrenocortical cells", "Chromaffin cells", "Endothelial cells",
                      "Fibroblasts", "Myeloid cells")
x = readRDS("data/sn.PPGLs.sct.latent_var_nonreg2.rds")
if(sum(x$orig.ident == sample & x$Cell_Supertype == "Sustentacular cells") != 1 ){
  x = x[,which((x$orig.ident == sample & x$Cell_Supertype %in% c("Tumour", "Sustentacular cells")) | (x$Cell_Supertype %in% normal.cell.types))]
} else {
  x = x[,which((x$orig.ident == sample & x$Cell_Supertype %in% c("Tumour")) | (x$Cell_Supertype %in% normal.cell.types))]
}
gene.position.file = "data/gene_positions.tsv"
gene.positions = read.delim(gene.position.file, sep="\t", header=F, as.is=T)
mat = as.matrix(x@assays$RNA@counts[match(gene.positions[,1], rownames(x@assays$RNA@counts)),])
annot = x@meta.data[,c("Cell_Supertype"), drop=F]
rm(x)
gc()

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=mat,
                                    annotations_file=annot,
                                    delim="\t",
                                    gene_order_file=gene.position.file,
                                    ref_group_names=normal.cell.types)
outdir = file.path("Results", sample, "InferCNV_updated")
if(!dir.exists(outdir)){dir.create(outdir)}
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=outdir,
                             cluster_by_groups=FALSE, 
                             denoise=TRUE,
                             HMM=TRUE,
                             analysis_mode = "subclusters",
                             tumor_subcluster_partition_method = "qnorm",
                             tumor_subcluster_pval = 0.01,
                             num_threads=16)
saveRDS(infercnv_obj, file=file.path(outdir, "output.rds"))

