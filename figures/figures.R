# Setup
source("colour_scheme.R")

# Dingbats mess with Adobe Illustrator apparently
pdf = function(...){grDevices::pdf(..., useDingbats = FALSE)}

# utility function for boxplots
add.counts = function(x,y=""){levels(x) = paste0(levels(x), " (n = ", table(x[!is.na(y)]), ")"); return(x)}

# ggplot setup
theme_set(theme_bw() + theme(panel.grid = element_blank())) # plain theme
# makes legend spacing smaller and removes legend box
t = theme(legend.background = element_blank(),
          legend.spacing = unit(0.1, "mm"),
          legend.key.size = unit(2, "mm"))
# Makes colour legend big and square
t2 = guides(colour = guide_legend(override.aes = list(shape = 15, size=6)))
t3 = guides(colour = guide_legend(override.aes = list(shape = 15, size=6, ncol=2)))
# theme for removing panel border
rmborder = list(scale_x_continuous(expand=expand_scale(mult = c(0.0, 0.2))), scale_y_continuous(expand=expand_scale(mult=c(0.1, 0.1))), theme(panel.border = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), legend.position = c(1, 0.5), legend.justification = "left", plot.margin = unit(c(1,1,1,1), "inches"), panel.spacing = unit(c(1, 1, 1, 1), "inches")), xlab(""), ylab(""), coord_fixed(ratio=1, clip="off"))
# honestly not sure I still need this
# supposed to make margins even between different plots with different sized legends
inclborder = list(scale_x_continuous(expand=expand_scale(mult = c(0.05, -(1-(1.05/1.2)))), limits = function(x){c(x[1], x[2]+0.2*(x[2]-x[1]))}),
                  #scale_y_continuous(expand=expand_scale(mult = rep(mean(c(0.05, -(1-(1.05/1.2)))), 2)), limits = function(x){c(x[1] - 0.1*(x[2]-x[1]), x[2]+0.1*(x[2]-x[1]))}),
                  theme(legend.position = c(1.05, 0.5), legend.justification = "left", plot.margin = unit(c(1, 3, 1, 1), "inches"), panel.spacing = unit(c(1, 1, 1, 1), "inches")))
# enforces square ggplot panel in real units
square.ratio = function(x){
  range2 = function(x){sum(c(-1,1)*range(x, na.rm=T))}
  return(coord_fixed(clip="off", ratio = range2(x[,1])/range2(x[,2])))
}




########################## Figure 1 #############################
fig1 = "Raw figures/Figure 1/Components"
dir.create(fig1, recursive=T, showWarnings = F)

#genotype.order = c("Normal", "RET", "HRAS", "NF1", "TMEM127", "MAX", "H3F3A", "MAML3", "FH", "EPAS1", "VHL", "SDHA", "SDHB", "SDHD", "SDHB_HN", "SDHD_HN", "Unknown") # old cluster labels

genotype.order = c("Normal",  "SDHA", "SDHB", "SDHC", "SDHD", "SETD2", "SDHB_HN", "SDHD_HN", "VHL", "EGLN1", "EPAS1", "FH", "RET", "HRAS", "NF1", "BRAF", "TMEM127", "NGFR", "MAX", "CSDE1", "H3F3A", "KIF1B", "IDH1", "MAML3", "Unknown")

sn.PPGLs$Genotype = factor(sn.PPGLs$Genotype %>% as.character, levels = genotype.order %>% intersect(unique(sn.PPGLs$Genotype)))
sn.PPGLs$Genotype2 = factor(sn.PPGLs$Genotype2 %>% as.character, levels = genotype.order %>% intersect(unique(sn.PPGLs$Genotype2)))
bulk.meta$Genotype = factor(as.character(bulk.meta$Genotype), levels = intersect(genotype.order, unique(bulk.meta$Genotype)))
levels(sn.PPGLs$Subgroup3)[levels(sn.PPGLs$Subgroup3) == "MAML"] = "MAML3"

#bulk.meta$OldCluster = bulk.meta$Cluster
bulk.meta$Cluster = all.clusters2[match(bulk.meta$OldCluster, all.clusters)]
bulk.meta$Cluster = factor(bulk.meta$Cluster, levels=c("Normal", "C1Ai", "C1Aii", "C1Bi", "C1Bii", "C2A", "C2Bi", "C2Bii", "C2C"))
sn.PPGLs$Subgroup4 = factor(all.clusters2[match(sn.PPGLs$Subgroup3, all.clusters)], levels=levels(bulk.meta$Cluster))

cols = intersect(colnames(pseudobulk.meta), colnames(sn.PPGLs@meta.data))
pseudobulk.meta[,cols] = sn.PPGLs@meta.data[match(pseudobulk.meta$orig.ident, sn.PPGLs$orig.ident),cols]
pseudobulk.meta$Subgroup4 = sn.PPGLs$Subgroup4[match(pseudobulk.meta$orig.ident, sn.PPGLs$orig.ident)]
pseudobulk.meta$Malignancy = c("Benign" = "Benign", "Metastatic" = "Malignant", "Benign(Local invasion)" = "Benign")[pseudobulk.meta$Malignancy_RAW]
pseudobulk.meta$Subgroup2[pseudobulk.meta$orig.ident=="A7"] = "IDC"
pseudobulk.meta$Subgroup3[pseudobulk.meta$orig.ident=="A7"] = "IDC"
pseudobulk.meta$Subgroup4[pseudobulk.meta$orig.ident=="A7"] = "C2Bi"
sn.umap.ratio = square.ratio(Embeddings(sn.PPGLs, reduction="umap"))

# SN genotype pie chart (Panel B)
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )
g = table(sn.PPGLs$Genotype[!duplicated(sn.PPGLs$Sample)])
pdf(file.path(fig1, "B_genotype_pie_chart.pdf"))
ggplot(data.frame(Genotype=factor(names(g), levels=levels(sn.PPGLs$Genotype)), Count=as.numeric(g)), aes(x=factor(""), fill=Genotype, y=Count))+
  geom_col(width = 1)+
  coord_polar("y") +  blank_theme +
  theme(axis.text.x=element_blank()) + scale_fill_manual(name="", values=genotype.cols, labels = paste0(names(g), " (n=", g, ")")) + guides(fill = guide_legend(override.aes = list(shape = 15, size=6), ncol=2))
dev.off()
# SN UMAP (Panel C)

pt.size=0.1
x=cbind(sn.PPGLs@meta.data, Embeddings(sn.PPGLs, reduction='umap'))
pdf(file.path(fig1, "C_sn_UMAP_genotype.pdf"))
ggplot() + geom_point(data=x, mapping=aes(x=UMAP_1, y=UMAP_2), col=ifelse(sn.PPGLs$Cell_Type %in% c("Tumour", "Chromaffin cells"), NA, alpha("grey80", 0.05)),size=pt.size, show.legend = FALSE) + geom_point(data=x, mapping=aes(x=UMAP_1, y=UMAP_2, col=factor(ifelse(Cell_Type %in% c("Tumour", "Chromaffin cells"), as.character(Genotype), NA), levels=levels(sn.PPGLs$Genotype))), size=pt.size) + sn.umap.ratio + scale_color_manual(name="", values=genotype.cols[levels(sn.PPGLs$Genotype)], labels=sapply(levels(sn.PPGLs$Genotype), function(x){if(x %in% c("Unknown","Normal")){x}else{bquote(italic(.(x)))}}), breaks=levels(sn.PPGLs$Genotype)) + t + guides(colour = guide_legend(override.aes = list(shape = 15, size=6, ncol=2)))
dev.off()
pdf(file.path(fig1, "C_sn_UMAP_cell_type.pdf"))
#ggplot() + geom_point(data=x, mapping=aes(x=UMAP_1, y=UMAP_2, col=Cell_Type), size=pt.size) + sn.umap.ratio +
#  scale_color_manual(name="", values=cell.cols[levels(sn.PPGLs$Cell_Type)], labels=replace(paste0(c("", paste0("(", substr(levels(sn.PPGLs$Cell_Type)[-1],1,1), ") ")), levels(sn.PPGLs$Cell_Type)), match("Macrophages/Monocytes", levels(sn.PPGLs$Cell_Type)), expression("(M"~Phi~") Macrophages/Monocytes")), breaks=levels(sn.PPGLs$Cell_Type))+
#  t + t3 + geom_text(data=x %>% filter(!Cell_Type %in% c("Macrophages/Monocytes", "Tumour")) %>% group_by(Cell_Type) %>% summarise(UMAP_1=median(UMAP_1), UMAP_2=median(UMAP_2), label=substr(Cell_Type[1],1,1)) %>% ungroup, mapping=aes(x=UMAP_1, y=UMAP_2, label=label)) + geom_text(data = x %>% filter(Cell_Type == "Macrophages/Monocytes") %>% group_by(Cell_Type) %>% summarise(UMAP_1=median(UMAP_1), UMAP_2 = median(UMAP_2)) %>% ungroup, mapping=aes(x=UMAP_1, y=UMAP_2), label=expression("M"~Phi), parse=T) + theme(legend.text.align=0)
ggplot() + geom_point(data=x, mapping=aes(x=UMAP_1, y=UMAP_2, col=Cell_Type), size=pt.size) + sn.umap.ratio +
  scale_color_manual(name="", values=cell.cols[levels(sn.PPGLs$Cell_Type)], labels=replace(paste0(c("", paste0("(", substr(levels(sn.PPGLs$Cell_Type)[-1],1,1), ") ")), levels(sn.PPGLs$Cell_Type)), match("Myeloid cells", levels(sn.PPGLs$Cell_Type)), expression("(My) Myeloid cells")), breaks=levels(sn.PPGLs$Cell_Type))+
  t + t3 + geom_text(data=x %>% filter(!Cell_Type %in% c("Myeloid cells", "Tumour")) %>% group_by(Cell_Type) %>% summarise(UMAP_1=median(UMAP_1), UMAP_2=median(UMAP_2), label=substr(Cell_Type[1],1,1)) %>% ungroup, mapping=aes(x=UMAP_1, y=UMAP_2, label=label)) + geom_text(data = x %>% filter(Cell_Type == "Myeloid cells") %>% group_by(Cell_Type) %>% summarise(UMAP_1=median(UMAP_1), UMAP_2 = median(UMAP_2)) %>% ungroup, mapping=aes(x=UMAP_1, y=UMAP_2), label=expression("My"), parse=T) + theme(legend.text.align=0)

dev.off()
#pdf(file.path(out.dir, "C_sn_UMAP_cell_type_labels_offset.pdf"))
#ggplot() + geom_point(data=x, mapping=aes(x=UMAP_1, y=UMAP_2, col=Cell_Type), size=pt.size) + sn.umap.ratio +
#  scale_color_manual(name="", values=cell.cols[levels(sn.PPGLs$Cell_Type)], labels=replace(paste0(c("", paste0("(", substr(levels(sn.PPGLs$Cell_Type)[-1],1,1), ") ")), levels(sn.PPGLs$Cell_Type)), match("Myeloid cells", levels(sn.PPGLs$Cell_Type)), expression("(My) Myeloid cells")), breaks=levels(sn.PPGLs$Cell_Type))+
#  t + t3 + geom_text(data=x %>% filter(!Cell_Type %in% c("Myeloid cells", "Tumour")) %>% group_by(Cell_Type) %>% summarise(UMAP_1=median(UMAP_1)-1, UMAP_2=median(UMAP_2)+1, label=substr(Cell_Type[1],1,1)) %>% ungroup, mapping=aes(x=UMAP_1, y=UMAP_2, label=label)) + geom_text(data = x %>% filter(Cell_Type == "Myeloid cells") %>% group_by(Cell_Type) %>% summarise(UMAP_1=median(UMAP_1)-1, UMAP_2 = median(UMAP_2)+1) %>% ungroup, mapping=aes(x=UMAP_1, y=UMAP_2), label=expression("My"), parse=T) + theme(legend.text.align=0)

#dev.off()
pdf(file.path(fig1, "C_sn_UMAP_cell_type_no_letters.pdf"))
#ggplot() + geom_point(data=x, mapping=aes(x=UMAP_1, y=UMAP_2, col=Cell_Type), size=pt.size) + sn.umap.ratio +
#  scale_color_manual(name="", values=cell.cols[levels(sn.PPGLs$Cell_Type)], labels=replace(paste0(c("", paste0("(", substr(levels(sn.PPGLs$Cell_Type)[-1],1,1), ") ")), levels(sn.PPGLs$Cell_Type)), match("Macrophages/Monocytes", levels(sn.PPGLs$Cell_Type)), expression("(M"~Phi~") Macrophages/Monocytes")), breaks=levels(sn.PPGLs$Cell_Type))+
#   t + t3  +theme(legend.text.align=0)
ggplot() + geom_point(data=x, mapping=aes(x=UMAP_1, y=UMAP_2, col=Cell_Type), size=pt.size) + sn.umap.ratio +
  scale_color_manual(name="", values=cell.cols[levels(sn.PPGLs$Cell_Type)], labels=replace(paste0(c("", paste0("(", substr(levels(sn.PPGLs$Cell_Type)[-1],1,1), ") ")), levels(sn.PPGLs$Cell_Type)), match("Myeloid cells", levels(sn.PPGLs$Cell_Type)), expression("(My) Myeloid cells")), breaks=levels(sn.PPGLs$Cell_Type))+
  t + t3  +theme(legend.text.align=0)

dev.off()
pdf(file.path(fig1, "C_sn_UMAP_sample_id.pdf"))
ggplot() + geom_point(data=x, mapping=aes(x=UMAP_1, y=UMAP_2, col=Sample), size=pt.size) + sn.umap.ratio + t + guides(colour = guide_legend(override.aes = list(shape = 15, size=6, ncol=2)))
dev.off()
# Dotplot (Panel D)

dotplot.genes = read.delim("Raw tables/Pseudobulk/Pseudobulk_one_vs_rest_key_genes_with_bulk_annotation - Pseudobulk_one_vs_rest_key_genes_with_bulk_annotation.tsv", header=T, as.is=T, sep='\t')

x = do.call(cbind, lapply(unique(sn.PPGLs$Sample), function(s){rowSums(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Sample==s & ((sn.PPGLs$Cell_Type %in% c("Chromaffin cells") & grepl("NAM", sn.PPGLs$orig.ident))|(sn.PPGLs$Cell_Type=="Tumour"))])})) %>%
  cbind(do.call(cbind, lapply(levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")), function(cell_type){
    rowSums(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Cell_Type==cell_type])
  }))) %>% LogNormalize(scale.factor=1e6) # logCPM
colnames(x) = c(unique(sn.PPGLs$Sample), levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")))
y = do.call(cbind, lapply(unique(sn.PPGLs$Sample), function(s){rowMeans(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Sample==s & ((sn.PPGLs$Cell_Type %in% c("Chromaffin cells") & grepl("NAM", sn.PPGLs$orig.ident))|(sn.PPGLs$Cell_Type=="Tumour"))] > 0)})) %>%
  cbind(do.call(cbind, lapply(levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")), function(cell_type){
    rowMeans(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Cell_Type==cell_type] > 0)
  })))# logCPM
colnames(y) = c(unique(sn.PPGLs$Sample), levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")))
z = data.frame(Sample=colnames(y), stringsAsFactors = F)
z$Cell_Type = ifelse(z$Sample %in% sn.PPGLs$Cell_Type, z$Sample, ifelse(z$Sample %in% sn.PPGLs$Sample[which(sn.PPGLs$Genotype=="Normal")], "Chromaffin cells", "Tumour"))
f = c(which(z$Cell_Type=="Tumour")[order(sn.PPGLs$Genotype2[match(z$Sample[which(z$Cell_Type=="Tumour")], sn.PPGLs$Sample)])],
      which(z$Cell_Type=="Chromaffin cells"),
      which(z$Cell_Type %in% setdiff(levels(sn.PPGLs$Cell_Type), c("Tumour", "Chromaffin cells")))[order(sn.PPGLs$Cell_Type[match(z$Sample[which(z$Cell_Type %in% setdiff(levels(sn.PPGLs$Cell_Type), c("Tumour", "Chromaffin cells")))], sn.PPGLs$Cell_Type)])])
z=z[f,]
x=x[,f]
y=y[,f]
z$Genotype = sn.PPGLs$Genotype[match(z$Sample, sn.PPGLs$Sample)]
z$Cluster = sn.PPGLs$Subgroup4[match(z$Sample, sn.PPGLs$Sample)]


#genes = c(unlist(lapply(setdiff(na.omit(unique(z$Cluster)), "Normal"), function(cluster){
#  cluster = gsub("SDHx\\ \\(H.*$", "SDHx_HN", cluster)
#  pseudobulk.tumour.vs.chromaffin.de.by.subtype %>% filter(!grepl("-|LINC|^RP", gene), Cluster==cluster, logFC > 0, adj.P.Val < 0.05, DE.same.direction.vs.rest.of.same.samples, DE.same.direction.vs.rest.of.tumours) %>% head(2) %>% pull(gene)
#})))
#genes = c(genes, unlist(lapply(na.omit(z$Cell_Type), function(cell_type){
#  if(cell_type == "Sustentacular cells"){
#    return(c("SOX10", "CDH19", "S100A1"))
#  } else {
#    return(head(pseudobulk.normal.signatures[[make.names(cell_type)]], 2))
#  }
#}))) %>% unique()

genes = dotplot.genes$gene[dotplot.genes$Figure=="1F"] %>% unique()

x=x[genes,]
y=y[genes,]
z[,"Cell Type"] = z$`Cell_Type`
hm.annot = HeatmapAnnotation(df = z[,c("Genotype", "Cluster", "Cell Type")], which="column", col = list("Genotype" = genotype.cols, "Cluster" = cluster.cols, "Cell Type" = cell.cols), show_legend = FALSE)
hm1 = HeatmapDotPlot(colour = x, size = y,
               scale = TRUE, cell.size = 0.5,
               cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
               col=colorRamp2(c(-3,0,3),c("blue","grey90","red")),
               name="Scaled expression", row_names_gp=gpar(fontface="italic"))

f2 = bulk.ccp[[9]]$consensusTree$order
f3 = f2[order(factor(as.character(bulk.meta$Cluster[f2]), levels=c(setdiff(levels(bulk.meta$Cluster), c("Cortical admixture", "Normal", "C2C")), "C2C", "Normal")))]


hm.bulk.annot = HeatmapAnnotation(df = bulk.meta[f3,c("Genotype", "Cluster")], col =list("Genotype" = genotype.cols, Cluster = cluster.cols), show_annotation_name = F)
hm2 = Heatmap(bulk.mat2[genes,f3] %>% t %>% scale %>% t, bottom_annotation=hm.bulk.annot, col = colorRamp2(c(-3,0,3),c("blue","grey90","red")), show_heatmap_legend = FALSE, cluster_columns = F, cluster_rows=FALSE, show_column_names = FALSE, show_row_names = F, row_names_gp = gpar(fontface="italic"))

#pdf(file.path(fig1, "F_dotplot_heatmaps_together.pdf"), width=18, height=9)
#hm2+hm1
#dev.off()

#pdf(file.path(fig1, "F_bulk_heatmap.pdf"))
#hm.bulk.annot = HeatmapAnnotation(df = bulk.meta[f3,c("Genotype", "Cluster")], col =list("Genotype" = genotype.cols, Cluster = cluster.cols), show_annotation_name = T)
#hm2 = Heatmap(bulk.mat2[genes,f3] %>% t %>% scale %>% t, name="Scaled expression", bottom_annotation=hm.bulk.annot, col = colorRamp2(c(-3,0,3),c("blue","white","red")), show_heatmap_legend = TRUE, cluster_columns = F, cluster_rows=FALSE, show_column_names = FALSE, show_row_names = T, row_names_gp = gpar(fontface="italic"))
#hm2
#dev.off()

pdf(file.path(fig1, "F_sn_dotplot.pdf"))
hm.annot = HeatmapAnnotation(df = z[,c("Genotype", "Cluster", "Cell Type")], which="column", col = list("Genotype" = genotype.cols, "Cluster" = cluster.cols, "Cell Type" = cell.cols), show_legend = TRUE)
hm1 = HeatmapDotPlot(colour = x, size = y, name="Scaled expression",
                     scale = TRUE, cell.size = 0.5,
                     cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
                     col=colorRamp2(c(-3,0,3),c("blue","grey90","red")), row_names_gp=gpar(fontface="italic"))
hm1
dev.off()

x = cbind(do.call(cbind,lapply(c("Adrenal pheochromocytoma", "AT-PGL", "HN-PGL")%>% setNames(c("Adrenal-PCC", "AT-PGL", "HN-PGL")), function(loc){
  rowSums(sn.PPGLs@assays$RNA@counts[,which(sn.PPGLs$Cell_Type=="Tumour" & sn.PPGLs$Location==loc)])
})),
do.call(cbind,lapply(levels(sn.PPGLs$Cell_Type) %>% setdiff("Tumour") %>% (function(x){setNames(x,x)}), function(cell_type){
  rowSums(sn.PPGLs@assays$RNA@counts[,which(sn.PPGLs$Cell_Type==cell_type)])
}))
) %>% as.matrix() %>% LogNormalize(scale.factor=1e6)
y = cbind(do.call(cbind,lapply(c("Adrenal pheochromocytoma", "AT-PGL", "HN-PGL")%>% setNames(c("Adrenal-PCC", "AT-PGL", "HN-PGL")), function(loc){
  rowMeans((sn.PPGLs@assays$RNA@counts[,which(sn.PPGLs$Cell_Type=="Tumour" & sn.PPGLs$Location==loc)] %>% as.matrix)> 0)
})),
do.call(cbind,lapply(levels(sn.PPGLs$Cell_Type) %>% setdiff("Tumour") %>% (function(x){setNames(x,x)}), function(cell_type){
  rowMeans((sn.PPGLs@assays$RNA@counts[,which(sn.PPGLs$Cell_Type==cell_type)] %>% as.matrix()) > 0)
}))) %>% as.matrix()

z = data.frame(`Cell Type` = colnames(x) %>% (function(x){factor(x,levels=unique(x))}),stringsAsFactors = F, check.names=F)
hm.annot = HeatmapAnnotation(df = z, which="column", col = list("Cell Type" = c(setNames(c(head(location.cols,3), cell.cols[colnames(x)[-(1:3)]]), colnames(x)))), show_legend = FALSE, show_annotation_name = F)

hm = HeatmapDotPlot(colour = x[genes,], size = y[genes,],
                     scale = TRUE, cell.size = 0.5,
                     cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
                     col=colorRamp2(c(-3,0,3),c("blue","grey90","red")),
                     name="Scaled expression", row_names_gp=gpar(fontface="italic"))

pdf(file.path(fig1, "F_dotplot_location_grouped.pdf"))
draw(hm)
dev.off()


genes = c('TH', genes)
locs = c("Adrenal pheochromocytoma", "AT-PGL", "HN-PGL", "Normal") %>% setNames(c("Adrenal-PCC", "AT-PGL", "HN-PGL", "Normal"))
cells = levels(sn.PPGLs$Cell_Type)[-1]
hmlist = lapply(locs, function(loc){
  x = do.call(cbind, lapply(cells, function(cell_type){
    rowSums(sn.PPGLs@assays$RNA@counts[,which(sn.PPGLs$Location == loc & gsub("Tumour", "Chromaffin cells", as.character(sn.PPGLs$Cell_Type)) == cell_type & (!(sn.PPGLs$Cell_Type=="Chromaffin cells" & loc != "Normal"))), drop=F])
  })
  ) %>% as.matrix() %>% LogNormalize(scale.factor=1e6)
  y = do.call(cbind,lapply(cells, function(cell_type){
    rowMeans((sn.PPGLs@assays$RNA@counts[,which(sn.PPGLs$Location == loc & gsub("Tumour", "Chromaffin cells", as.character(sn.PPGLs$Cell_Type)) == cell_type & (!(sn.PPGLs$Cell_Type=="Chromaffin cells" & loc != "Normal"))), drop=F] %>% as.matrix) > 0)
  }))
  colnames(x) = cells
  colnames(y) = cells
  z = data.frame(`Cell Type` = cells %>% (function(x){factor(x,levels=unique(x))}),stringsAsFactors = F, check.names=F)
  hm.annot = HeatmapAnnotation(df = z, which="column", col = list("Cell Type" = cell.cols), show_legend = FALSE, show_annotation_name = F)
  if(loc != "Normal"){
    colnames(x)[1] = names(locs)[locs==loc]
  } else {
    colnames(x)[1] = "Chromaffin cells"
  }
  y[is.na(y)] = 0
  hm = HeatmapDotPlot(colour = x[genes,] %>% as.matrix, size = y[genes,] %>% as.matrix,
                      scale = FALSE, cell.size = 0.5,
                      cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
                      col=colorRamp2(c(0, 4, 8),c("blue","grey90","red")),
                      name="Expression", row_names_gp=gpar(fontface="italic"),
                      show_column_names = loc==tail(locs,1), show_row_names = TRUE,
                      show_heatmap_legend=loc==locs[1],
                      row_title = names(locs)[locs==loc], row_title_side="left", row_title_rot = 0)
  return(hm)
})
gc(full=TRUE)
source("htlistfix.R")
hm = Reduce((function(x,y){x %v% y}), hmlist)

pdf(file.path(fig1, "F_dotplot_location_grouped_normals_split.pdf"), height=10)
draw(hm)
dev.off()

# Cell type heatmap (Panel E)
y = factor(sn.PPGLs$Cell_Type,
           levels = c("Tumour", setdiff(levels(sn.PPGLs$Cell_Type), "Tumour")))
y2 = sn.PPGLs$Sample

z = table(y, y2)
y = data.frame(Sample=sort(unique(sn.PPGLs$Sample)), stringsAsFactors = F)
y$Genotype = factor(sn.PPGLs$Genotype[match(y$Sample, sn.PPGLs$Sample)],
                    levels = c("Normal",
                               "RET",
                               "HRAS",
                               "NF1",
                               "TMEM127",
                               "MAML3",
                               "MAX",
                               "H3F3A",
                               "FH",
                               "EPAS1",
                               "VHL",
                               "SDHA",
                               "SDHB",
                               "SDHD",
                               "Unknown"))
y$Sample = factor(y$Sample, levels = y$Sample[order(y$Genotype)])
y = y[order(y$Genotype),]

w = do.call('cbind', lapply(levels(y$Sample), function(i){z[,i]}))
colnames(w) = levels(y$Sample)
w = 1000*sweep(w, 2, colSums(w), "/")
w = as.data.frame(w)
w = do.call('cbind', lapply(w, round, digits=0))
rownames(w) = rownames(z)
for(i in 1:ncol(w)){
  j = w[,i]
  w[j==max(j),i] = w[j==max(j),i] - (sum(j) - 1000)
}
z = data.frame(do.call('cbind', lapply(as.data.frame(w), function(x){do.call('c', lapply(1:nrow(w), function(i){rep(rownames(w)[i], x[i])}))})))
z = as.matrix(z)


pdf(file.path(fig1, "E_cell_type_distribution_barplot_with_normals.pdf"),width=10)
annot = HeatmapAnnotation(df = y[,"Genotype",drop=F], name = "Genotype", col = list("Genotype"=genotype.cols), show_annotation_name = F, show_legend = F)
hm = Heatmap(z[,], name="Cell Type", col=cell.cols, bottom_annotation = annot, heatmap_legend_param = list(title="Cell Type", ncol=1), show_heatmap_legend = T)
draw(hm, annotation_legend_list = Legend(labels = levels(y$Genotype), title = "Genotype", type = "grid", 
                                        labels_gp = gpar(fontface = ifelse(levels(y$Genotype) %in% c("Normal", "Unknown"), "plain", "italic")),
                                        legend_gp = gpar(fill = genotype.cols[levels(y$Genotype)])))
dev.off()
pdf(file.path(fig1, "E_cell_type_distribution_barplot_with_normals_raster.pdf"),width=10)
annot = HeatmapAnnotation(df = y[,"Genotype",drop=F], name = "Genotype", col = list("Genotype"=genotype.cols), show_annotation_name = F, show_legend = F)
hm = Heatmap(z[,], name="Cell Type", col=cell.cols, bottom_annotation = annot, heatmap_legend_param = list(title="Cell Type", ncol=1), show_heatmap_legend = T, use_raster = TRUE)
draw(hm, annotation_legend_list = Legend(labels = levels(y$Genotype), title = "Genotype", type = "grid", 
                                         labels_gp = gpar(fontface = ifelse(levels(y$Genotype) %in% c("Normal", "Unknown"), "plain", "italic")),
                                         legend_gp = gpar(fill = genotype.cols[levels(y$Genotype)])))
dev.off()
f=which(!colnames(z) %in% c("E240", "E243"))
pdf(file.path(fig1, "E_cell_type_distribution_barplot_no_normals.pdf"),width=10)
annot = HeatmapAnnotation(df = y[f,"Genotype",drop=F], name = "Genotype", col = list("Genotype"=genotype.cols), show_annotation_name = F, show_legend = F)
hm = Heatmap(z[,f], name="Cell Type", col=cell.cols, bottom_annotation = annot, heatmap_legend_param = list(title="Cell Type", ncol=1), show_heatmap_legend = T)
draw(hm, annotation_legend_list = Legend(labels = levels(y$Genotype)[-1], title = "Genotype", type = "grid", 
                                         labels_gp = gpar(fontface = ifelse(levels(y$Genotype)[-1] %in% c("Normal", "Unknown"), "plain", "italic")),
                                         legend_gp = gpar(fill = genotype.cols[levels(y$Genotype)[-1]])))
dev.off()
pdf(file.path(fig1, "E_cell_type_distribution_barplot_no_normals_raster.pdf"),width=10)
annot = HeatmapAnnotation(df = y[f,"Genotype",drop=F], name = "Genotype", col = list("Genotype"=genotype.cols), show_annotation_name = F, show_legend = F)
hm = Heatmap(z[,f], name="Cell Type", col=cell.cols, bottom_annotation = annot, heatmap_legend_param = list(title="Cell Type", ncol=1), show_heatmap_legend = T, use_raster = TRUE)
draw(hm, annotation_legend_list = Legend(labels = levels(y$Genotype)[-1], title = "Genotype", type = "grid", 
                                         labels_gp = gpar(fontface = ifelse(levels(y$Genotype)[-1] %in% c("Normal", "Unknown"), "plain", "italic")),
                                         legend_gp = gpar(fill = genotype.cols[levels(y$Genotype)[-1]])))
dev.off()

# Panel F - InferCNV copy number
source("htlistfix.R")
library(RColorBrewer)
s = c("VHP35T" = "data/PPGL_SNP6/AF_VPH35T_(CytoScanHD_Array)",
      "VHP09T" = "data/PPGL_SNP6/V-PH-09T_(CytoScanHD_Array)/",
      "VPH20T" = "data/PPGL_SNP6/V-PH-20T_(CytoScanHD_Array)/",
      "PGL3" = "data/PPGL_SNP6/10_119_P018_P1-1_(CytoScanHD_Array)",
      "PGL1" = "data/PPGL_SNP6/10_119_P018_P1-1_(CytoScanHD_Array)",
      "A7" = NA, "E210" = NA)
#names(s) = c("07AH391", "VCBP14T", "A11", "PGL3", "PGL1")
infercnv_objs = lapply(names(s), function(s) readRDS(file.path("Results", s, "InferCNV_updated/output.rds"))) %>% setNames(names(s))

#infercnv_obj = infercnv_objs[[1]]
infercnv_sust_expr = lapply(system(intern=T, command = "ls Results/*/InferCNV_updated/output.rds"), function(f){
  infercnv_obj = readRDS(f)
  y = t(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_observations$labels] - 1)
  y = y[sn.PPGLs@meta.data[rownames(y),"Cell_Supertype"] %in% "Sustentacular cells",]
  return(y)
})
genes = Reduce(intersect, lapply(infercnv_sust_expr, colnames))
infercnv_sust_expr = do.call(rbind, lapply(infercnv_sust_expr, function(x){x[,genes]}))
gc(full=TRUE)
# InferCNV subclusters
sample="PGL1"
infercnv_obj = infercnv_objs[[sample]]
genes = rownames(infercnv_obj@expr.data) %>% intersect(colnames(infercnv_sust_expr))
genes.chr = infercnv_obj@gene_order[genes,]$chr %>% (function(x){factor(x,levels=unique(x))})

x = t(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_references$labels][,infercnv_obj@tumor_subclusters$hc$all_references$order] -1)
x = rbind(x[,genes], infercnv_sust_expr[,genes])
x = x[order(sn.PPGLs@meta.data[rownames(x), "Cell_Supertype"]),]
y = t(infercnv_obj@expr.data[genes,infercnv_obj@tumor_subclusters$hc$all_observations$labels] - 1)
y = y[sn.PPGLs@meta.data[rownames(y),"Cell_Supertype"] %in% "Tumour",]

x = x[,genes]
y = y[,genes]

#set.seed(1); fx=sort(sample(which(sn.PPGLs@meta.data[rownames(x), "Cell_Type"] != "Sustentacular cells"), 1000))
set.seed(1); fx=sort(sample(1:nrow(x), 1000))
col_fun1 = colorRamp2(breaks = c(-.1,0,.1), colors=c("blue", "white","red"))
col_fun2 = colorRamp2(c(-0.3,0,0.3), c("blue","white","red"))
x2 = HeatmapAnnotation(df = data.frame("Cell Type" = as.character(sn.PPGLs@meta.data[rownames(x)[fx], "Cell_Supertype"]), check.names = F), col=list("Cell Type" = cell.cols), which="row", show_legend = T, show_annotation_name = TRUE, annotation_name_side = "top")
hm1 = Heatmap(x[fx,],left_annotation = x2, use_raster = TRUE, cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=F, col = col_fun1, show_heatmap_legend = F)

#rawcopy.segs = read.delim("data/PPGL_SNP6//segments.txt", header=T, as.is=T, sep='\t')
rawcopy.segs = read.delim(file.path(s[sample], "segments.txt"), header=T, as.is=T, sep='\t')

row.order = colnames(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_observations$labels])[infercnv_obj@tumor_subclusters$hc$all_observations$order] %>% intersect(rownames(y))

y2 = HeatmapAnnotation(df = data.frame(`Tumour Cluster` = as.character(cutree(infercnv_obj@tumor_subclusters$hc$all_observations, k=2)), check.names=F), col=list("Tumour Cluster" = setNames(brewer.pal(3, "RdYlGn")[c(1,3)], c("1", "2"))),which="row", show_annotation_name = TRUE, annotation_name_side = "bottom", show_legend=TRUE)
#y2 = HeatmapAnnotation(df = data.frame(`Tumour Cluster` = as.character(cutree(infercnv_obj@tumor_subclusters$hc$all_observations, k=2)[match(rownames(y), infercnv_obj@tumor_subclusters$hc$all_observations$labels)][match(row.order, rownames(y))]), check.names=F), col=list("Tumour Cluster" = setNames(brewer.pal(3, "RdYlGn")[c(1,3)], c("1", "2"))),which="row", show_annotation_name = TRUE, annotation_name_side = "bottom", show_legend=TRUE)

y3 = HeatmapAnnotation(which="column", df = data.frame(SNP6 = sapply(genes, function(i){
  a=infercnv_obj@gene_order[i,,drop=F]
  b=rawcopy.segs[which(rawcopy.segs$Chromosome==a[,1]),,drop=F]
  return(tail(b$Value[b$Start < a$stop], 1))
}), Chromosome=infercnv_obj@gene_order[genes,]$chr)[,c("SNP6", "Chromosome"), drop=F], col = list(SNP6 = col_fun2, Chromosome = setNames(rep(brewer.pal(n=3, "Accent")[2:3], 11), unique(infercnv_obj@gene_order$chr))), show_annotation_name = T, annotation_name_side = "right", show_legend = c(TRUE, FALSE))
hm2=Heatmap(y, name="InferCNV", left_annotation=y2, bottom_annotation=y3, show_heatmap_legend = TRUE, col=col_fun1, cluster_columns = FALSE, cluster_rows = infercnv_obj@tumor_subclusters$hc$all_observations, show_row_names = F, show_column_names=F, use_raster=TRUE)
#hm2=Heatmap(y[row.order,], name="InferCNV", left_annotation=y2, bottom_annotation=y3, show_heatmap_legend = TRUE, col=col_fun1, cluster_columns = FALSE, cluster_rows = F, show_row_names = F, show_column_names=F, use_raster=TRUE)

pdf(file.path(fig1, "F_infercnv_heatmap_with_sust_PGL1.pdf"))
draw(hm1 %v% hm2)
dev.off()

sample="VHP35T"
infercnv_obj = infercnv_objs[[sample]]
genes = rownames(infercnv_obj@expr.data) %>% intersect(colnames(infercnv_sust_expr))
genes.chr = infercnv_obj@gene_order[genes,]$chr %>% (function(x){factor(x,levels=unique(x))})

x = t(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_references$labels][,infercnv_obj@tumor_subclusters$hc$all_references$order] -1)
x = rbind(x[,genes], infercnv_sust_expr[,genes])
x = x[order(sn.PPGLs@meta.data[rownames(x), "Cell_Supertype"]),]
y = t(infercnv_obj@expr.data[genes,infercnv_obj@tumor_subclusters$hc$all_observations$labels] - 1)
y = y[sn.PPGLs@meta.data[rownames(y),"Cell_Supertype"] %in% "Tumour",]

x = x[,genes]
y = y[,genes]

#set.seed(1); fx=sort(sample(which(sn.PPGLs@meta.data[rownames(x), "Cell_Type"] != "Sustentacular cells"), 1000))
set.seed(1); fx=sort(sample(1:nrow(x), 1000))
col_fun1 = colorRamp2(breaks = c(-.1,0,.1), colors=c("blue", "white","red"))
col_fun2 = colorRamp2(c(-0.3,0,0.3), c("blue","white","red"))
x2 = HeatmapAnnotation(df = data.frame("Cell Type" = as.character(sn.PPGLs@meta.data[rownames(x)[fx], "Cell_Supertype"]), check.names = F), col=list("Cell Type" = cell.cols), which="row", show_legend = T, show_annotation_name = TRUE, annotation_name_side = "top")
hm1 = Heatmap(x[fx,],left_annotation = x2, use_raster = TRUE, cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=F, col = col_fun1, show_heatmap_legend = F)

#rawcopy.segs = read.delim("data/PPGL_SNP6//segments.txt", header=T, as.is=T, sep='\t')
rawcopy.segs = read.delim(file.path(s[sample], "segments.txt"), header=T, as.is=T, sep='\t')

row.order = colnames(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_observations$labels])[infercnv_obj@tumor_subclusters$hc$all_observations$order] %>% intersect(rownames(y))

#y2 = HeatmapAnnotation(df = data.frame(`Tumour Cluster` = as.character(cutree(infercnv_obj@tumor_subclusters$hc$all_observations, k=2)), check.names=F), col=list("Tumour Cluster" = setNames(brewer.pal(3, "RdYlGn")[c(1,3)], c("1", "2"))),which="row", show_annotation_name = TRUE, annotation_name_side = "bottom", show_legend=TRUE)
#y2 = HeatmapAnnotation(df = data.frame(`Tumour Cluster` = as.character(cutree(infercnv_obj@tumor_subclusters$hc$all_observations, k=2)[match(rownames(y), infercnv_obj@tumor_subclusters$hc$all_observations$labels)][match(row.order, rownames(y))]), check.names=F), col=list("Tumour Cluster" = setNames(brewer.pal(3, "RdYlGn")[c(1,3)], c("1", "2"))),which="row", show_annotation_name = TRUE, annotation_name_side = "bottom", show_legend=TRUE)

y3 = HeatmapAnnotation(which="column", df = data.frame(SNP6 = sapply(genes, function(i){
  a=infercnv_obj@gene_order[i,,drop=F]
  b=rawcopy.segs[which(rawcopy.segs$Chromosome==a[,1]),,drop=F]
  return(tail(b$Value[b$Start < a$stop], 1))
}), Chromosome=infercnv_obj@gene_order[genes,]$chr)[,c("SNP6", "Chromosome"), drop=F], col = list(SNP6 = col_fun2, Chromosome = setNames(rep(brewer.pal(n=3, "Accent")[2:3], 11), unique(infercnv_obj@gene_order$chr))), show_annotation_name = T, annotation_name_side = "right", show_legend = c(TRUE, FALSE))
#hm2=Heatmap(y, name="InferCNV", bottom_annotation=y3, show_heatmap_legend = TRUE, col=col_fun1, cluster_columns = FALSE, cluster_rows = infercnv_obj@tumor_subclusters$hc$all_observations, show_row_names = F, show_column_names=F, use_raster=TRUE)
hm2=Heatmap(y[row.order,], name="InferCNV",  bottom_annotation=y3, show_heatmap_legend = TRUE, col=col_fun1, cluster_columns = FALSE, cluster_rows = F, show_row_names = F, show_column_names=F, use_raster=TRUE)

pdf(file.path(fig1, "F_infercnv_heatmap_with_sust_VPH35T.pdf"))
draw(hm1 %v% hm2)
dev.off()

for(sample in unique(sn.PPGLs$orig.ident)){
  print(sample)
  infercnv_file = file.path("Results", sample, "InferCNV_updated/output.rds")
  if(file.exists(infercnv_file)){
  infercnv_obj = readRDS(infercnv_file)
  genes = rownames(infercnv_obj@expr.data) %>% intersect(colnames(infercnv_sust_expr))
  genes.chr = infercnv_obj@gene_order[genes,]$chr %>% (function(x){factor(x,levels=unique(x))})
  
  x = t(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_references$labels][,infercnv_obj@tumor_subclusters$hc$all_references$order] -1)
  x = rbind(x[,genes], infercnv_sust_expr[,genes])
  x = x[order(sn.PPGLs@meta.data[rownames(x), "Cell_Supertype"]),]
  y = t(infercnv_obj@expr.data[genes,infercnv_obj@tumor_subclusters$hc$all_observations$labels] - 1)
  y = y[sn.PPGLs@meta.data[rownames(y),"Cell_Supertype"] %in% "Tumour",]
  
  x = x[,genes]
  y = y[,genes]
  
  #set.seed(1); fx=sort(sample(which(sn.PPGLs@meta.data[rownames(x), "Cell_Type"] != "Sustentacular cells"), 1000))
  set.seed(1); fx=sort(sample(1:nrow(x), 1000))
  col_fun1 = colorRamp2(breaks = c(-.1,0,.1), colors=c("blue", "white","red"))
  col_fun2 = colorRamp2(c(-0.3,0,0.3), c("blue","white","red"))
  x2 = HeatmapAnnotation(df = data.frame("Cell Type" = as.character(sn.PPGLs@meta.data[rownames(x)[fx], "Cell_Supertype"]), check.names = F), col=list("Cell Type" = cell.cols), which="row", show_legend = T, show_annotation_name = TRUE, annotation_name_side = "top")
  hm1 = Heatmap(x[fx,],left_annotation = x2, use_raster = TRUE, cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=F, col = col_fun1, show_heatmap_legend = F)
  
  #rawcopy.segs = read.delim("data/PPGL_SNP6//segments.txt", header=T, as.is=T, sep='\t')
  # don't use
  rawcopy.segs = read.delim(file.path(s[1], "segments.txt"), header=T, as.is=T, sep='\t')
  
  row.order = colnames(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_observations$labels])[infercnv_obj@tumor_subclusters$hc$all_observations$order] %>% intersect(rownames(y))
  
  #y2 = HeatmapAnnotation(df = data.frame(`Tumour Cluster` = as.character(cutree(infercnv_obj@tumor_subclusters$hc$all_observations, k=2)), check.names=F), col=list("Tumour Cluster" = setNames(brewer.pal(3, "RdYlGn")[c(1,3)], c("1", "2"))),which="row", show_annotation_name = TRUE, annotation_name_side = "bottom", show_legend=TRUE)
  #y2 = HeatmapAnnotation(df = data.frame(`Tumour Cluster` = as.character(cutree(infercnv_obj@tumor_subclusters$hc$all_observations, k=2)[match(rownames(y), infercnv_obj@tumor_subclusters$hc$all_observations$labels)][match(row.order, rownames(y))]), check.names=F), col=list("Tumour Cluster" = setNames(brewer.pal(3, "RdYlGn")[c(1,3)], c("1", "2"))),which="row", show_annotation_name = TRUE, annotation_name_side = "bottom", show_legend=TRUE)
  
  y3 = HeatmapAnnotation(which="column", df = data.frame(SNP6 = sapply(genes, function(i){
    a=infercnv_obj@gene_order[i,,drop=F]
    b=rawcopy.segs[which(rawcopy.segs$Chromosome==a[,1]),,drop=F]
    return(tail(b$Value[b$Start < a$stop], 1))
  }), Chromosome=infercnv_obj@gene_order[genes,]$chr)[,c("Chromosome"), drop=F], col = list(SNP6 = col_fun2, Chromosome = setNames(rep(brewer.pal(n=3, "Accent")[2:3], 11), unique(infercnv_obj@gene_order$chr))), show_annotation_name = T, annotation_name_side = "right", show_legend = c(FALSE))
  #hm2=Heatmap(y, name="InferCNV", bottom_annotation=y3, show_heatmap_legend = TRUE, col=col_fun1, cluster_columns = FALSE, cluster_rows = infercnv_obj@tumor_subclusters$hc$all_observations, show_row_names = F, show_column_names=F, use_raster=TRUE)
  hm2=Heatmap(y[row.order,], name="InferCNV",  bottom_annotation=y3, show_heatmap_legend = TRUE, col=col_fun1, cluster_columns = FALSE, cluster_rows = F, show_row_names = F, show_column_names=F, use_raster=TRUE)
  
  pdf(file.path(fig1, paste0("F_infercnv_heatmap_with_sust_", sample, ".pdf")))
  draw(hm1 %v% hm2)
  dev.off()
  }
}

# 4 in one
s = s[c(1:3,5)]
infercnv_objs = infercnv_objs[c(1:3,5)]

genes = Reduce(intersect, lapply(infercnv_objs, function(x){rownames(x@expr.data)})) %>% intersect(colnames(infercnv_sust_expr))
genes.chr = infercnv_objs[[1]]@gene_order[genes,]$chr %>% (function(x){factor(x,levels=unique(x))})

hm.list = lapply(names(s), function(sample){
#print(sample)
infercnv_obj = infercnv_objs[[sample]]
x = t(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_references$labels][,infercnv_obj@tumor_subclusters$hc$all_references$order] -1)
x = x[order(sn.PPGLs@meta.data[rownames(x), "Cell_Supertype"]),]


#y = do.call(rbind, lapply(infercnv_objs, function(infercnv_obj) t(infercnv_obj@expr.data[genes,infercnv_obj@tumor_subclusters$hc$all_observations$labels] - 1)))
y = t(infercnv_obj@expr.data[genes,infercnv_obj@tumor_subclusters$hc$all_observations$labels] - 1)
#f = match(colnames(sn.PPGLs)[which(sn.PPGLs$orig.ident == sample & sn.PPGLs$Cell_Type=="Sustentacular cells")], rownames(x))
#y = rbind(y[,], x[f,genes])
#x=x[-f,genes]
x = x[,genes]
y = y[,genes]

#set.seed(1); fx=sort(sample(which(sn.PPGLs@meta.data[rownames(x), "Cell_Type"] != "Sustentacular cells"), 1000))
set.seed(1); fx=sort(sample(1:nrow(x), 1000))
col_fun1 = colorRamp2(breaks = c(-.1,0,.1), colors=c("blue", "white","red"))
col_fun2 = colorRamp2(c(-0.3,0,0.3), c("blue","white","red"))
x2 = HeatmapAnnotation(df = data.frame("Cell Type" = sn.PPGLs@meta.data[rownames(x)[fx], "Cell_Supertype"], check.names = F), col=list("Cell Type" = cell.cols), which="row", show_legend = F, show_annotation_name = TRUE, annotation_name_side = "top")
hm1 = Heatmap(x[fx,],left_annotation = x2, use_raster = TRUE, cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=F, col = col_fun1, show_heatmap_legend = F)

#rawcopy.segs = read.delim("data/PPGL_SNP6//segments.txt", header=T, as.is=T, sep='\t')
rawcopy.segs = read.delim(file.path(s[sample], "segments.txt"), header=T, as.is=T, sep='\t')
#y2 = HeatmapAnnotation(df = data.frame(`Tumour Cluster` = as.character(cutree(infercnv_obj@tumor_subclusters$hc$all_observations, k=2)), check.names=F), col=list("Tumour Cluster" = setNames(brewer.pal(3, "RdYlGn")[c(1,3)], c("1", "2"))),which="row", show_annotation_name = TRUE, annotation_name_side = "bottom", show_legend=TRUE)
#y2 = HeatmapAnnotation(df = data.frame(`Cell Type` = sn.PPGLs@meta.data[rownames(y),"Cell_Supertype"], check.names=F), col=list("Cell Type" = cell.cols[c("Tumour", "Sustentacular cells", as.character(unique(sn.PPGLs@meta.data[rownames(x)[fx],"Cell_Supertype"])))]),which="row", show_annotation_name = FALSE, annotation_name_side = "bottom", show_legend=TRUE)

cols.keep = c("SNP6")
if(sample == tail(names(s), 1)){cols.keep = c("SNP6", "Chromosome")}
#cols.keep = "Chromosome"
y3 = HeatmapAnnotation(which="column", df = data.frame(SNP6 = sapply(genes, function(i){
  a=infercnv_obj@gene_order[i,,drop=F]
  b=rawcopy.segs[which(rawcopy.segs$Chromosome==a[,1]),,drop=F]
  return(tail(b$Value[b$Start < a$stop], 1))
}), Chromosome=infercnv_obj@gene_order[genes,]$chr)[,cols.keep, drop=F], col = list(SNP6 = col_fun2, Chromosome = setNames(rep(brewer.pal(n=3, "Accent")[2:3], 11), unique(infercnv_obj@gene_order$chr))), show_annotation_name = T, annotation_name_side = "left", show_legend = c(sample==tail(names(s),1), F)[1:(length(cols.keep))])
#hm2=Heatmap(y, name="InferCNV", left_annotation=y2, bottom_annotation=y3, show_heatmap_legend = TRUE, col=col_fun1, cluster_columns = FALSE, cluster_rows = infercnv_obj@tumor_subclusters$hc$all_observations, show_row_names = F, show_column_names=F, use_raster=TRUE)
if(sample == tail(names(s),1) | "SNP6" %in% cols.keep){
  hm2=Heatmap(y, name="InferCNV", row_title=sn.PPGLs$Sample[match(sample, sn.PPGLs$orig.ident)], bottom_annotation=y3, show_heatmap_legend = TRUE, col=col_fun1, cluster_columns = FALSE, cluster_rows = T, show_row_dend = F, show_row_names = F, show_column_names=F, use_raster=TRUE, column_split = genes.chr, column_title_gp = gpar(col=NA), row_title_rot = 0)
} else {
  hm2=Heatmap(y, name="InferCNV", row_title=sn.PPGLs$Sample[match(sample, sn.PPGLs$orig.ident)], show_heatmap_legend = TRUE, col=col_fun1, cluster_columns = FALSE, cluster_rows = T, show_row_dend = F, show_row_names = F, show_column_names=F, use_raster=TRUE, column_split = genes.chr, column_title_gp = gpar(col=NA), row_title_rot = 0)
}

return(hm2)
})
sample=names(s)[1]
infercnv_obj = infercnv_objs[[1]]
x = t(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_references$labels][,infercnv_obj@tumor_subclusters$hc$all_references$order] -1)
x = rbind(x[,genes], infercnv_sust_expr[,genes])
x = x[order(sn.PPGLs@meta.data[rownames(x), "Cell_Supertype"]),]


#y = do.call(rbind, lapply(infercnv_objs, function(infercnv_obj) t(infercnv_obj@expr.data[genes,infercnv_obj@tumor_subclusters$hc$all_observations$labels] - 1)))
y = t(infercnv_obj@expr.data[genes,infercnv_obj@tumor_subclusters$hc$all_observations$labels] - 1)
#f = match(colnames(sn.PPGLs)[which(sn.PPGLs$orig.ident == sample & sn.PPGLs$Cell_Type=="Sustentacular cells")], rownames(x))
#y = rbind(y, x[f,genes])
#x=x[-f,genes]
y = y[,genes]
x = x[,genes]

#set.seed(1); fx=sort(sample(which(sn.PPGLs@meta.data[rownames(x), "Cell_Supertype"] != "Sustentacular cells"), 1000))
set.seed(1); fx=sort(sample(1:nrow(x), 4000))
col_fun1 = colorRamp2(breaks = c(-.1,0,.1), colors=c("blue", "white","red"))
col_fun2 = colorRamp2(c(-0.3,0,0.3), c("blue","white","red"))
x2 = HeatmapAnnotation(df = data.frame("Cell Type" = sn.PPGLs@meta.data[rownames(x)[fx], "Cell_Supertype"], check.names = F), col=list("Cell Type" = cell.cols), which="row", show_legend = F, show_annotation_name = TRUE, annotation_name_side = "top")
hm1 = Heatmap(x[fx,], left_annotation = x2, use_raster = TRUE, cluster_rows=F, cluster_columns=F, show_row_names=F, show_column_names=F, col = col_fun1, show_heatmap_legend = F, column_split = genes.chr, column_title_gp = gpar(col=NA))

hm_list = Reduce(function(x,y){x %v% y}, c(list(hm1), hm.list))

l = intersect(levels(sn.PPGLs$Cell_Supertype), c("Sustentacular cells", names(infercnv_objs[[1]]@reference_grouped_cell_indices)))
lgd_list = list(
  Legend(labels=l,
         legend_gp = gpar(fill=cell.cols[l]),
         background = "#FFFFFF", title = "Cell Type")
)
#hm_list = hm1 %v% hm2

pdf(file.path(fig1, "F_infercnv_heatmap_4_samples_with_SNP6.pdf"))
draw(hm_list, heatmap_legend_list = lgd_list)
dev.off()

infercnv_obj = infercnv_objs$PGL3
x = t(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_references$labels][,infercnv_obj@tumor_subclusters$hc$all_references$order] -1)
x = rbind(x, t(infercnv_obj@expr.data[,infercnv_obj@tumor_subclusters$hc$all_observations$labels][,infercnv_obj@tumor_subclusters$hc$all_observations$order] -1))

x = x[order(sn.PPGLs@meta.data[rownames(x), "Cell_Supertype"]),]
x = x[sn.PPGLs@meta.data[rownames(x), "orig.ident"] == "PGL3",]
genes = colnames(x)
genes.chr = infercnv_obj@gene_order[genes,]$chr %>% (function(x){factor(x,levels=unique(x))})

rawcopy.segs = read.delim(file.path("data/PPGL_SNP6/10_119_P018_P1-1_(CytoScanHD_Array)", "segments.txt"), header=T, as.is=T, sep='\t')

y3 = HeatmapAnnotation(which="column", df = data.frame(SNP6 = sapply(genes, function(i){
  a=infercnv_obj@gene_order[i,,drop=F]
  b=rawcopy.segs[which(rawcopy.segs$Chromosome==a[,1]),,drop=F]
  return(tail(b$Value[b$Start < a$stop], 1))
}), Chromosome=infercnv_obj@gene_order[genes,]$chr)[,c("SNP6", "Chromosome"), drop=F],
col = list(SNP6 = col_fun2,
           Chromosome = setNames(rep(brewer.pal(n=3, "Accent")[2:3], 11),
                                 unique(infercnv_obj@gene_order$chr))),
show_annotation_name = T, annotation_name_side = "left", show_legend = c(FALSE))
y4 = HeatmapAnnotation(which="column", df = data.frame(SNP6 = sapply(genes, function(i){
  a=infercnv_obj@gene_order[i,,drop=F]
  b=rawcopy.segs[which(rawcopy.segs$Chromosome==a[,1]),,drop=F]
  return(tail(b$Value[b$Start < a$stop], 1))
}), Chromosome=infercnv_obj@gene_order[genes,]$chr)[,c("Chromosome"), drop=F],
col = list(SNP6 = col_fun2,
           Chromosome = setNames(rep(brewer.pal(n=3, "Accent")[2:3], 11),
                                 unique(infercnv_obj@gene_order$chr))),
show_annotation_name = T, annotation_name_side = "left", show_legend = c(FALSE))
x2 = HeatmapAnnotation(df = data.frame("Cell Type" = sn.PPGLs@meta.data[rownames(x), "Cell_Supertype"], check.names = F),
                       col=list("Cell Type" = cell.cols), which="row", show_legend = T, show_annotation_name = F, annotation_name_side = "top")
hm1 = Heatmap(x, left_annotation = x2, bottom_annotation=y3, use_raster = TRUE, cluster_rows=T, cluster_columns=F, show_row_names=F, show_column_names=F, show_row_dend=F, row_title_rot=0, col = col_fun1, show_heatmap_legend = F, column_split = genes.chr, column_title_gp = gpar(col=NA), row_split = sn.PPGLs@meta.data[rownames(x),"Cell_Supertype"])
hm2 = Heatmap(x, left_annotation = x2, bottom_annotation=y4, use_raster = TRUE, cluster_rows=T, cluster_columns=F, show_row_names=F, show_column_names=F, show_row_dend=F, row_title_rot=0, col = col_fun1, show_heatmap_legend = F, column_split = genes.chr, column_title_gp = gpar(col=NA), row_split = sn.PPGLs@meta.data[rownames(x),"Cell_Supertype"])

pdf(file.path(fig1, "F_infercnv_single_sample.pdf"))
print(hm1)
print(hm2)
dev.off()

########################## Figure 2 #############################
fig2 = "Raw figures/Figure 2/Components"
bulk.umap.ratio = square.ratio(bulk.meta[,c("UMAP_1", "UMAP_2")])
dir.create(fig2, showWarnings = F, recursive = TRUE)
size=4
pdf(file.path(fig2, "bulk_UMAP_genotype.pdf"), width=10, height=7)
ggplot(bulk.meta %>% arrange(Genotype != "Unknown"), aes(x=UMAP_1, y=UMAP_2, col=Genotype)) + geom_point(size=size) + genotype.scale + t + t2 + bulk.umap.ratio
dev.off()

pdf(file.path(fig2, "bulk_UMAP_old_clusters.pdf"), width=10, height=7)
ggplot(bulk.meta, aes(x=UMAP_1, y=UMAP_2, col=OldCluster)) + geom_point(size=size) + cluster.scale + bulk.umap.ratio + t + t2
dev.off()

pdf(file.path(fig2, "bulk_UMAP_clusters.pdf"), width=10, height=7)
ggplot(bulk.meta, aes(x=UMAP_1, y=UMAP_2, col=Cluster)) + geom_point(size=size) + cluster.scale + bulk.umap.ratio + t + t2
dev.off()

pdf(file.path(fig2, "bulk_UMAP_anatomic_location.pdf"), width=10, height=7)
ggplot(bulk.meta %>% arrange(!is.na(Location)) %>% mutate(Location=ifelse(is.na(Location), "Unknown", Location)), aes(x=UMAP_1, y=UMAP_2, col=Location)) + geom_point(size=size) + bulk.umap.ratio+ t + t2 + scale_colour_manual(values = location.cols)
dev.off()

genotype.table = table(bulk.meta$Genotype, bulk.meta$Sample)
genotype.table = as.data.frame(do.call(cbind, lapply(bulk.meta$Sample, function(s){genotype.table[,s]})))
colnames(genotype.table) = bulk.meta$Sample
top.annot = HeatmapAnnotation(df = bulk.meta[,c("Cluster","Location"), drop=F], col = list("Cluster" = cluster.cols, "Location"= location.cols, "Malignancy" = malignancy.cols), show_annotation_name = T)
hm = Heatmap(genotype.table %>% as.matrix, col = c("gray90", "black"), show_column_names = F, show_row_names=T, column_split = NULL, cluster_column_slices = F, cluster_columns = bulk.ccp[[9]]$consensusTree, cluster_rows = FALSE, top_annotation = top.annot, show_heatmap_legend = FALSE, border = FALSE)

pdf(file.path(fig2, "genotype_heatmap_consensus_clustering.pdf"))
print(hm)
dev.off()

f = order(bulk.meta$Genotype)
f2 = order(bulk.meta$Cluster[f])
top.annot = HeatmapAnnotation(df = bulk.meta[f,][f2,][,c("Cluster","Location"), drop=F], col = list("Cluster" = cluster.cols, "Location"= location.cols, "Malignancy" = malignancy.cols), show_annotation_name = T)
hm = Heatmap(genotype.table[,f][,f2] %>% as.matrix, col = c("gray90", "black"), show_column_names = F, show_row_names=T, column_split = bulk.meta$Cluster[f][f2], cluster_column_slices = F, cluster_columns = F, cluster_rows = FALSE, top_annotation = top.annot, show_heatmap_legend = FALSE, border = FALSE)

pdf(file.path(fig2, "genotype_heatmap_genotype_ordered.pdf"), width=12)
print(hm)
dev.off()

x = data.frame(Sample = c(bulk.meta$Sample, pseudobulk.meta$Sample),
               Genotype = c(as.character(bulk.meta$Genotype), as.character(pseudobulk.meta$Genotype)) %>% factor(levels = genotype.order %>% intersect(c(bulk.meta$Genotype %>% as.character, pseudobulk.meta$Genotype %>% as.character))),
               Cluster = c(bulk.meta$Cluster %>% as.character, pseudobulk.meta$Subgroup4 %>% as.character) %>% factor(levels=levels(bulk.meta$Cluster)),
               Location=c(bulk.meta$Location, c("Adrenal pheochromocytoma" = "Adrenal", "AT-PGL" = "Extraadrenal", "HN-PGL" = "Head and neck", "Normal" = NA)[pseudobulk.meta$Location]),
               Malignancy = c(bulk.meta$Malignancy, pseudobulk.meta$Malignancy), stringsAsFactors = F)
y = table(x$Genotype, x$Sample)
y = do.call(cbind, lapply(1:ncol(y), function(i){z = data.frame(a=y[,i]); colnames(z) = colnames(y)[i];return(z)}))
for(i in unique(pseudobulk.meta$Sample)){y[,i]=y[,i]*2}

y = y[,x$Sample]
f = order(x$Genotype)
f2 = order(x$Cluster[f])
top.annot = HeatmapAnnotation(df = x[f,][f2,][,c("Cluster","Location"), drop=F], col = list("Cluster" = cluster.cols, "Location"= location.cols, "Malignancy" = malignancy.cols), show_annotation_name = T)
hm = Heatmap(y[,f][,f2] %>% as.matrix, col = c("gray90", "black", "red"), show_column_names = F, show_row_names=T, column_split = x$Cluster[f][f2], cluster_column_slices = F, cluster_columns = F, cluster_rows = FALSE, top_annotation = top.annot, show_heatmap_legend = FALSE, border = FALSE, row_names_gp = gpar(fontface = ifelse(rownames(y[,f][,f2]) %in% c("Normal", "Unknown"), "plain", "italic")))
lgd_list = list(
  Legend(labels=c("Bulk", "Single-nuclei"),
         legend_gp = gpar(fill=c("black", "red")),
         background = "#FFFFFF", title = "Type")
)

pdf(file.path(fig2, "genotype_heatmap_genotype_ordered_with_pseudobulk.pdf"), width=12)
draw(hm, heatmap_legend_list = lgd_list,merge_legends=TRUE)
dev.off()

x = data.frame(Sample = c(bulk.meta$Sample, pseudobulk.meta$Sample),
               Genotype = c(as.character(bulk.meta$Genotype), as.character(pseudobulk.meta$Genotype)) %>% factor(levels = genotype.order %>% intersect(c(bulk.meta$Genotype %>% as.character, pseudobulk.meta$Genotype %>% as.character))),
               Cluster = c(bulk.meta$Cluster %>% as.character, pseudobulk.meta$Subgroup4 %>% as.character) %>% factor(levels=levels(bulk.meta$Cluster)),
               Location=c(bulk.meta$Location, c("Adrenal pheochromocytoma" = "Adrenal", "AT-PGL" = "Extraadrenal", "HN-PGL" = "Head and neck", "Normal" = NA)[pseudobulk.meta$Location]),
               Malignancy = c(bulk.meta$Malignancy, pseudobulk.meta$Malignancy), stringsAsFactors = F)
y = table(x$Genotype, x$Sample)
y = do.call(cbind, lapply(1:ncol(y), function(i){z = data.frame(a=y[,i]); colnames(z) = colnames(y)[i];return(z)}))
for(i in unique(pseudobulk.meta$Sample)){y[,i]=y[,i]*2}

y = y[,x$Sample]
f = order(x$Genotype)
f2 = order(x$Cluster[f])
left.annot = HeatmapAnnotation(df = x[f,][f2,][,c("Cluster","Location"), drop=F], col = list("Cluster" = cluster.cols, "Location"= location.cols, "Malignancy" = malignancy.cols), show_annotation_name = T, which='row')
hm = Heatmap(y[,f][,f2] %>% as.matrix %>% t, col = c("gray90", "black", "red"), show_column_names = T, show_row_names=F, row_split = x$Cluster[f][f2], cluster_row_slices = F, cluster_rows = F, cluster_columns = FALSE, left_annotation = left.annot, show_heatmap_legend = FALSE, border = FALSE, row_title_rot = 0, column_names_gp = gpar(fontface = ifelse(rownames(y[,f][,f2]) %in% c("Normal", "Unknown"), "plain", "italic")))
lgd_list = list(
  Legend(labels=c("Bulk", "Single-nuclei"),
         legend_gp = gpar(fill=c("black", "red")),
         background = "#FFFFFF", title = "Type")
)

pdf(file.path(fig2, "genotype_heatmap_genotype_ordered_with_pseudobulk_transposed.pdf"))
draw(hm, heatmap_legend_list = lgd_list,merge_legends=TRUE)
dev.off()



pdf(file.path(fig2, "bulk_umap_with_pseudobulk_genotype.pdf"), width=10)
ggplot() + geom_point(data=bulk.meta %>% arrange(Genotype!="Unknown"), mapping=aes(x=UMAP_1, y=UMAP_2, col=Genotype), alpha=0.5, shape=16, size=size) + genotype.scale + bulk.umap.ratio + t + t2 + 
  geom_point(data = pseudobulk.meta, mapping = aes(x=Pseudobulk_UMAP_1, y=Pseudobulk_UMAP_2, fill=Genotype), size=size, shape=22) + guides(alpha=F, shape=F, fill=F) + genotype.scale2 + guides(col=guide_legend(override.aes=list(shape=15,size=6, alpha=c(rep(1,22),0.25))))
dev.off()

pdf(file.path(fig2, "bulk_umap_with_pseudobulk_genotype_no_alpha_change.pdf"), width=10)
ggplot() + geom_point(data=bulk.meta %>% arrange(Genotype!="Unknown"), mapping=aes(x=UMAP_1, y=UMAP_2, col=Genotype), alpha=1, shape=16, size=size) + genotype.scale + bulk.umap.ratio + t + t2 + 
  geom_point(data = pseudobulk.meta, mapping = aes(x=Pseudobulk_UMAP_1, y=Pseudobulk_UMAP_2, fill=Genotype), size=size, shape=22) + guides(alpha=F, shape=F, fill=F) + genotype.scale2 + guides(col=guide_legend(override.aes=list(shape=15,size=6, alpha=c(rep(1,22),0.25))))
dev.off()

pdf(file.path(fig2, "bulk_umap_pseudobulk_position_changes.pdf"), width=10)
ggplot() + geom_point(data=bulk.meta, mapping=aes(x=UMAP_1, y=UMAP_2), col=alpha("grey50", 0.1), size=size) +
  geom_point(data = pseudobulk.meta, mapping = aes(x=Pseudobulk_chromaffin_UMAP_1, y=Pseudobulk_chromaffin_UMAP_2, col="Pseudobulk (chromaffin-only)"), size=size) +
  geom_point(data = pseudobulk.meta, mapping = aes(x=Pseudobulk_UMAP_1, y=Pseudobulk_UMAP_2,  col="Pseudobulk (all cells)"), size=size) +
  geom_point(data = pseudobulk.meta, mapping = aes(x=Bulk_UMAP_1, y=Bulk_UMAP_2, col="Bulk"), size=size) +
  geom_segment(data=pseudobulk.meta, mapping = aes(x=Bulk_UMAP_1, xend=Pseudobulk_UMAP_1, y=Bulk_UMAP_2, yend=Pseudobulk_UMAP_2)) +
  geom_segment(data=pseudobulk.meta, mapping = aes(x=Pseudobulk_UMAP_1, xend=Pseudobulk_chromaffin_UMAP_1, y=Pseudobulk_UMAP_2, yend=Pseudobulk_chromaffin_UMAP_2)) +
  scale_colour_discrete(name="") + bulk.umap.ratio + t + t2
dev.off()

pdf(file.path(fig2, "bulk_UMAP_tcga_subtypes.pdf"), width=10, height=7)
ggplot(bulk.meta %>% arrange(!is.na(TCGA_Cluster)), aes(x=UMAP_1, y=UMAP_2, col=TCGA_Cluster)) + geom_point(size=size) + bulk.umap.ratio + t + t2 + scale_colour_discrete(name="TCGA Subtype")
dev.off()

pdf(file.path(fig2, "bulk_UMAP_comete_subtypes.pdf"), width=10, height=7)
ggplot(bulk.meta %>% arrange(!is.na(Castro_Vega_mRNA_Classification)), aes(x=UMAP_1, y=UMAP_2, col=Castro_Vega_mRNA_Classification)) + geom_point(size=size) + bulk.umap.ratio + t + t2 + scale_colour_discrete(name="COMETE Subtype")
dev.off()

pdf(file.path(fig2, "sn_UMAP_cell_types.pdf"))
for(reduction in c("Myeloid.cells", "Fibroblasts", "Endothelial.cells", "T.cells", "B.cells")){
  x = cbind(Embeddings(sn.PPGLs, reduction=reduction),sn.PPGLs@meta.data)
  colnames(x)[1:2] = c("UMAP_1", "UMAP_2")
  x = x %>% filter(!is.na(UMAP_1)) %>% mutate(Cell_Type = gsub("pos", "+", gsub("neg", "-", Cell_Type)))
  p = ggplot(x, aes(x=UMAP_1, y=UMAP_2, col=Cell_Type, shape=ifelse(Sample %in% c("E240", "E243"), "Normal", "Tumor"))) + geom_point(size=3) + square.ratio(x[,1:2]) + scale_color_manual(values=brewer.pal(9, "Set1")[-6], name="Cell type") + scale_shape_manual(name="", values=c("Tumor" = 16, "Normal" = 17))
  print(p + t + t2)
}
for(reduction in c("Fibroblasts", "Endothelial.cells")){
  x = cbind(Embeddings(sn.PPGLs, reduction=reduction),sn.PPGLs@meta.data)
  colnames(x)[1:2] = c("UMAP_1", "UMAP_2")
  x = x %>% filter(!is.na(UMAP_1)) %>% mutate(Cell_Type = gsub("pos", "+", gsub("neg", "-", Cell_Type)), Tumour = ifelse(Sample %in% c("E240", "E243"), "Normal", "Tumor"))
  y = x %>% group_by(Tumour) %>% summarise(mid_x=median(UMAP_1), mid_y = median(UMAP_2),
                                           mad_x = mad(UMAP_1), mad_y=mad(UMAP_2))
  p = ggplot(x, aes(col=Tumour, x=UMAP_1,y=UMAP_2)) + stat_ellipse(show.legend = F) + geom_point(size=3, mapping=aes(x=UMAP_1, y=UMAP_2, col=Cell_Type, shape=Tumour)) + square.ratio(x[,1:2]) + scale_color_manual(values=brewer.pal(9, "Set1")[-6] %>% setNames(c("Tumor", "Normal", sort(unique(x$Cell_Type)))), name="Cell type", breaks=sort(unique(x$Cell_Type))) + scale_shape_manual(name="", values=c("Tumor" = 16, "Normal" = 17))
  print(p + t + t2)
  p = ggplot(x, aes(col=Tumour, x=UMAP_1,y=UMAP_2)) + stat_ellipse(show.legend = F) + geom_point(size=3, mapping=aes(x=UMAP_1, y=UMAP_2, col=Cell_Type)) + square.ratio(x[,1:2]) + scale_color_manual(values=brewer.pal(9, "Set1")[-6] %>% setNames(c("Tumor", "Normal", sort(unique(x$Cell_Type)))), name="Cell type", breaks=sort(unique(x$Cell_Type))) + scale_shape_manual(name="", values=c("Tumor" = 16, "Normal" = 17))
  print(p + t + t2)
}
dev.off()

reduction="TplusB"
sn.PPGLs = sn.PPGLs %>% RunUMAPSubset(cells = which(sn.PPGLs$Cell_Supertype %in% c("B cells", "T/NK cells")), name=reduction, dims = 1:20, metric = "cosine", k = 15, resolution = 0.8)
pdf(file.path(fig2, "sn_UMAP_cell_types_immune_cells_together.pdf"))
x = cbind(Embeddings(sn.PPGLs, reduction=reduction),sn.PPGLs@meta.data)
colnames(x)[1:2] = c("UMAP_1", "UMAP_2")
x = x %>% filter(!is.na(UMAP_1)) %>% mutate(Cell_Type = gsub("pos", "+", gsub("neg", "-", Cell_Type)))
x$Cell_Type = factor(x$Cell_Type, levels=x %>% arrange(Cell_Supertype, Cell_Type) %>% pull(Cell_Type) %>% unique)
p = ggplot(x, aes(x=UMAP_1, y=UMAP_2, col=Cell_Type)) + geom_point(size=3) + square.ratio(x[,1:2]) + scale_color_manual(values=brewer.pal(9, "Set1")[-6], name="Cell type") #+ scale_shape_manual(name="", values=c("Tumor" = 16, "Normal" = 17))
print(p + t + t2)
dev.off()
########################## Figure 3 #############################
#A. GSVA enrichment of gene sets from the bulk for each cell type (cell type gene sets) from the SN. 
#B. Boxplots of enrichments of cell types in groups and immune fractions. Talking about interferons and PD1 expression. (Maybe replace boxplots with a ridge plot). 
#C. Leukocyte score by subtype from TCGA methylation. 
#D. CD274 RPPA from TCGA and mRNA from everything. (Brady bunch of ridge plots)
fig3 = "Raw figures/Figure 3/Components"
dir.create(fig3, recursive = T, showWarnings = F)
# A
bulk.gsva = read.delim("Raw tables/Bulk/bulk_gsva.tsv", sep='\t', header=T, row.names=1, check.names=F)

#mat =  as.matrix(bulk.gsva)[!grepl("HALLMARK", rownames(bulk.gsva)),]
#f = order(factor(rownames(mat), levels = c("Chromaffin cells", "Fibroblasts", "Endothelial cells", "Sustentacular cells", "Adrenocortical cells", "Macrophages/Monocytes", "T cells", "B cells", "Mast cells") %>% make.names()))
#rownames(mat) = gsub("\\.", " ", gsub("Macrophages.", "Macrophages/", rownames(mat)))
#f = order(factor(rownames(mat), levels = c("Chromaffin cells", "Fibroblasts", "Endothelial cells", "Sustentacular cells", "Adrenocortical cells", "Myeloid cells", "T/NK cells", "B cells", "Mast cells") %>% make.names()))
#rownames(mat) = gsub("\\.", " ", gsub("^T.", "T/", rownames(mat)))

mat = bulk.gsva[names(pseudobulk.normal.subtype.signatures),]
mat.split = sn.PPGLs$Cell_Supertype[match(rownames(mat), sn.PPGLs$Cell_Type %>% make.names)]
mat = mat[order(mat.split),]
mat.split = sort(mat.split)
cell.annot = HeatmapAnnotation(which='row', df=data.frame(`Cell Type` = mat.split, check.names=F), col=list("Cell Type" = cell.cols), show_legend = F, show_annotation_name = F)
rownames(mat) = gsub("\\.", " ", gsub("pos", "+", gsub("neg", "-", rownames(mat))))

pdf(file.path(fig3, "A_bulk_gsva_cell_types_with_chromaffin_pseudobulk_signatures.pdf"), width=10)
hm.annot = HeatmapAnnotation(which='column', df = bulk.meta[,"Cluster", drop=F], col = list("Cluster" = cluster.cols))
hm = Heatmap(mat, name = "GSVA score", top_annotation = hm.annot, left_annotation=cell.annot, show_column_names=F, row_split=mat.split, row_title_rot = 0)
hm
dev.off()
pdf(file.path(fig3, "A_bulk_gsva_cell_types_with_chromaffin_pseudobulk_signatures_consensus_cluster_order.pdf"), width=10)
hm.annot = HeatmapAnnotation(which='column', df = bulk.meta[,"Cluster", drop=F], col = list("Cluster" = cluster.cols))
hm = Heatmap(mat, name = "GSVA score", top_annotation = hm.annot, left_annotation=cell.annot, show_column_names=F, row_split=mat.split, row_title_rot=0, cluster_columns = bulk.ccp[[9]]$consensusTree)
hm
dev.off()
pdf(file.path(fig3, "A_bulk_gsva_cell_types_with_chromaffin_pseudobulk_signatures_consensus_split_by_cluster.pdf"), width=10)
hm.annot = HeatmapAnnotation(which='column', df = bulk.meta[,"Cluster", drop=F], col = list("Cluster" = cluster.cols), show_legend=F)
hm = Heatmap(mat, name = "GSVA score", top_annotation = hm.annot, left_annotation=cell.annot, show_column_names=F,row_split=mat.split, row_title_rot=0, cluster_columns = T, cluster_column_slices=F, column_title_rot=90, column_split = bulk.meta$Cluster, show_parent_dend_line = F)
hm
dev.off()

f2 = which(bulk.meta$Cluster %in% c("Normal", "Cortical admixture", "C2C"))
pdf(file.path(fig3, "A_bulk_gsva_cell_types_with_chromaffin_pseudobulk_signatures_no_normals.pdf"), width=10)
hm.annot = HeatmapAnnotation(which='column', df = bulk.meta[-f2,"Cluster", drop=F], col = list("Cluster" = cluster.cols))
hm = Heatmap(mat[,-f2], name = "GSVA score", top_annotation = hm.annot, show_column_names=F)
hm
dev.off()
pdf(file.path(fig3, "A_bulk_gsva_cell_types_with_chromaffin_pseudobulk_signatures_consensus_cluster_order_no_normals.pdf"), width=10)
f3 = setdiff(bulk.ccp[[9]]$consensusTree$order, f2)
hm.annot = HeatmapAnnotation(which='column', df = bulk.meta[f3,"Cluster", drop=F], col = list("Cluster" = cluster.cols))
hm = Heatmap(mat[,f3], name = "GSVA score", top_annotation = hm.annot, show_column_names=F, cluster_columns = F)
hm
dev.off()
pdf(file.path(fig3, "A_bulk_gsva_cell_types_with_chromaffin_pseudobulk_signatures_consensus_split_by_cluster_no_normals.pdf"), width=10)
hm.annot = HeatmapAnnotation(which='column', df = bulk.meta[-f2,"Cluster", drop=F], col = list("Cluster" = cluster.cols))
hm = Heatmap(mat[,-f2], name = "GSVA score", top_annotation = hm.annot, show_column_names=F, cluster_columns = T, column_split = sapply(as.integer(bulk.meta$Cluster[-f2]), function(i){paste(rep(" ", i), collapse="")}), show_parent_dend_line = F)
hm
dev.off()

#genes = unlist(lapply(unique(bulk.de.by.subtype$Cluster) %>% setdiff("Cortical.admixture"), function(x){
#  bulk.de.by.subtype %>% filter(logFC > 0, adj.P.Val <= 0.05, Cluster==x, gene %in% rownames(sn.PPGLs)) %>% pull(gene) %>% head(20)
#})) %>% unique()
g = pseudobulk.tumour.vs.chromaffin.de.by.subtype %>% group_by(Cluster) %>%
  arrange(Cluster) %>% 
  filter(DE.same.direction.vs.rest.of.same.samples) %>%
  filter(adj.P.Val <= 0.05) %>% arrange(-abs(logFC)) %>%
  mutate(Gene = gene, Sign = sign(logFC)) %>%
  filter(Gene %in% (bulk.de.by.subtype[bulk.de.by.subtype$Cluster==Cluster[1],] %>% 
                      filter(adj.P.Val <= 0.05) %>%
                      filter(sign(logFC) == Sign[match(gene, Gene)]) %>%
                      pull(gene))
                    ) %>%
  mutate(bulk.logFC = bulk.de.by.subtype$logFC[match(paste(gene, Cluster), (bulk.de.by.subtype %>% mutate(x=paste(gene,Cluster)) %>% pull(x)))]) %>%
  arrange(-abs(bulk.logFC)) %>%
  filter(gene %in% c(head(gene[logFC > 0], 10), head(gene[logFC < 0], 10))) %>% ungroup() %>% arrange(Cluster)
genes = g$gene %>% unique
g = g %>% filter(!duplicated(gene)) %>% cbind(table(g$gene, g$Cluster)[genes,] %>% (function(x){do.call(cbind, lapply(colnames(x) %>% setNames(colnames(x)), function(y){ifelse(x[,y] ==0, 0, ifelse(g$logFC[g$Cluster==y][match(rownames(x), g$gene[g$Cluster==y])] >0, "1", "-1"))}))}))

hm.annot = HeatmapAnnotation(which='column', df = bulk.meta[,"Cluster", drop=F], col = list("Cluster" = cluster.cols), show_legend=F, show_annotation_name=F)
#hm.annot2 = HeatmapAnnotation(which='row', df = data.frame(`Cell Signature` = sapply(genes, function(gene){x=names(pseudobulk.normal.signatures)[which(sapply(names(pseudobulk.normal.signatures), function(x){gene %in% (pseudobulk.normal.de %>% filter(Cluster==x, logFC >= 2, adj.P.Val <= 0.05) %>% pull(gene))}))]; if(length(x)==0){NA}else{gsub("\\.", " ", gsub("pos", "+", gsub("neg", "-", x[1])))}}), check.names=F), col = list("Cell Signature" = setNames(c30, gsub("\\.", " ", gsub("pos", "+", gsub("neg", "-", names(pseudobulk.normal.signatures))))) %>% (function(x){x[!is.na(names(x))]})), show_annotation_name = F)
hm.annot2 = HeatmapAnnotation(which='row', df=g[,(unique(pseudobulk.tumour.vs.chromaffin.de.by.subtype$Cluster) %>% sort)], col = (lapply(unique(pseudobulk.tumour.vs.chromaffin.de.by.subtype$Cluster) %>% (function(x) setNames(x,x)), function(x){c("-1" = "blue", "0" = "white", "1" = "red")})), show_legend = F)
hm = Heatmap(bulk.mat2[genes,] %>% t %>% scale %>% t, name="Scaled Expression", top_annotation=hm.annot,
             show_column_names=F, show_row_names=F, show_column_dend = F, 
             cluster_columns=T, cluster_column_slices = F,
             right_annotation=hm.annot2,
             column_split = bulk.meta$Cluster, show_parent_dend_line=F,
             #row_split = unlist(lapply(1:7, function(i){rep(paste0(rep(" ", i), collapse=''), 20)})),
             row_split = g$Cluster,
             row_title = " ", column_title_rot=90,
             cluster_row_slices = F)
#hm2 = Heatmap(pseudobulk.voom$E[genes,] %>% (function(x){sweep(x, 2, apply(x,2,min), "-")}), name = "Pseudobulk expression", column_split = gsub("pos", "+", gsub("neg", "-", gsub("\\.", " ", gsub(".*_", "", colnames(pseudobulk.voom$E))))), column_gap = unit(0,"cm"),
#              show_column_names=F, show_row_names=F, show_column_dend = F, cluster_rows = F, cluster_columns=T, cluster_column_slices=F, column_title_rot=90)
pdf(file.path(fig3, "bulk_top20_genes_heatmap_with_cluster_highlights.pdf"), width=20)
#hm + hm2
hm
dev.off()

# B boxplots

pdf(file.path(fig3, "B_bulk_gsva_violin.pdf"))
for(cell_type in setdiff(levels(sn.PPGLs$Cell_Type), "Tumour")){
  p = ggplot(bulk.meta %>% cbind(t(bulk.gsva)) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y=paste0("`", cell_type %>% make.names, "`"), fill="Cluster")) + geom_violin(show.legend = F) + geom_point(alpha=0.5, position=position_jitter(width=0.1), show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)] %>% unname)
  print(p + xlab("") + ylab(paste(gsub("s$", "", cell_type), "GSVA score")) + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
}
dev.off()

pdf(file.path(fig3, "B_bulk_gsva_violin_no_normals.pdf"))
for(cell_type in setdiff(levels(sn.PPGLs$Cell_Type), "Tumour")){
  p = ggplot(bulk.meta %>% cbind(t(bulk.gsva)) %>% filter(!Cluster %in% c("Normal", "Cortical admixture", "C2C")) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y=paste0("`", cell_type %>% make.names, "`"), fill="Cluster")) + geom_violin(show.legend = F) + geom_point(alpha=0.5, position=position_jitter(width=0.1), show.legend = F) + scale_fill_manual(values=cluster.cols[setdiff(levels(bulk.meta$Cluster), c("Normal", "Cortical admixture"))] %>% unname)
  print(p + xlab("") + ylab(paste(gsub("s$", "", cell_type), "GSVA score")) + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
}
dev.off()

pdf(file.path(fig3, "B_bulk_gsva_boxplot.pdf"))
for(cell_type in setdiff(levels(sn.PPGLs$Cell_Type), "Tumour")){
  p = ggplot(bulk.meta %>% cbind(t(bulk.gsva)) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y=paste0("`", cell_type %>% make.names, "`"), fill="Cluster")) + geom_boxplot(show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)] %>% unname)
  print(p + xlab("") + ylab(paste(gsub("s$", "", cell_type), "GSVA score")) + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
}
dev.off()

pdf(file.path(fig3, "B_bulk_gsva_boxplot_no_normals.pdf"))
for(cell_type in setdiff(levels(sn.PPGLs$Cell_Type), "Tumour")){
  p = ggplot(bulk.meta %>% cbind(t(bulk.gsva)) %>% filter(!Cluster %in% c("Normal", "Cortical admixture", "C2C")) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y=paste0("`", cell_type %>% make.names, "`"), fill="Cluster")) + geom_boxplot(show.legend = F) + scale_fill_manual(values=cluster.cols[setdiff(levels(bulk.meta$Cluster), c("Normal", "Cortical admixture"))] %>% unname)
  print(p + xlab("") + ylab(paste(gsub("s$", "", cell_type), "GSVA score")) + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
}
dev.off()

# C 
tcga.leuk = read.delim("Raw tables/Bulk/TCGA_all_leuk_estimate.masked.20170107.tsv", sep='\t', header=F, as.is=T) %>% filter(V1=="PCPG")
tcga.rppa = read.delim("Raw tables/Bulk/TCGA-PCPG-L4_RPPA.csv", sep=',', header=T, as.is=T)

bulk.meta$Leukocyte_score = tcga.leuk[match(bulk.meta$Sample, make.names(substr(tcga.leuk[,2], 1, 15))),3]

pdf(file.path(fig3, "C_tcga_leukocyte_score_boxplot.pdf"))
p = ggplot(bulk.meta %>% filter(!is.na(Leukocyte_score)) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="Leukocyte_score", fill="Cluster")) + geom_boxplot(show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)] %>% unname)
print(p + xlab("") + ylab("Leukocyte score") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
dev.off()
pdf(file.path(fig3, "C_tcga_leukocyte_score_boxplot_no_normals.pdf"))
p = ggplot(bulk.meta %>% filter(!is.na(Leukocyte_score))%>% filter(!Cluster %in% c("Normal", "Cortical admixture", "C2C")) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="Leukocyte_score", fill="Cluster")) + geom_boxplot(show.legend = F) + scale_fill_manual(values=cluster.cols[setdiff(levels(bulk.meta$Cluster), c("Normal", "Cortical admixture"))] %>% unname)
print(p + xlab("") + ylab("Leukocyte score") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
dev.off()

pdf(file.path(fig3, "C_tcga_leukocyte_score_violin.pdf"))
p = ggplot(bulk.meta %>% filter(!is.na(Leukocyte_score)) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="Leukocyte_score", fill="Cluster")) + geom_violin(show.legend = F) + geom_point(alpha=0.5, position=position_jitter(width = 0.1), show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)] %>% unname)
print(p + xlab("") + ylab("Leukocyte score") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
dev.off()
pdf(file.path(fig3, "C_tcga_leukocyte_score_violin_no_normals.pdf"))
p = ggplot(bulk.meta %>% filter(!is.na(Leukocyte_score))%>% filter(!Cluster %in% c("Normal", "Cortical admixture", "C2C")) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="Leukocyte_score", fill="Cluster")) + geom_violin(show.legend = F) + geom_point(alpha=0.5, position=position_jitter(width = 0.1), show.legend = F) + scale_fill_manual(values=cluster.cols[setdiff(levels(bulk.meta$Cluster), c("Normal", "Cortical admixture"))] %>% unname)
print(p + xlab("") + ylab("Leukocyte score") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
dev.off()

# D RPPA
f=match(bulk.meta$Sample, substr(tcga.rppa$Sample_ID,1,15) %>% make.names)
pdf(file.path(fig3, "D_tcga_pdl1_rppa_and_mrna_boxplot.pdf"))
p = ggplot(bulk.meta %>% mutate(PDL1=tcga.rppa[f,"PDL1"]) %>% filter(!is.na(PDL1)) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="PDL1", fill="Cluster")) + geom_boxplot(show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)[-1]] %>% unname)
print(p + xlab("") + ylab("CD274 RPPA (TCGA)") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
p = ggplot(bulk.meta %>% mutate(PDL1=bulk.mat2["CD274",] %>% as.numeric) %>% filter(!is.na(PDL1)) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="PDL1", fill="Cluster")) + geom_boxplot(show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)] %>% unname)
print(p + xlab("") + ylab("CD274 expression") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
dev.off()
pdf(file.path(fig3, "D_tcga_pdl1_rppa_and_mrna_boxplot_no_normals.pdf"))
p = ggplot(bulk.meta %>% mutate(PDL1=tcga.rppa[f,"PDL1"]) %>% filter(!is.na(PDL1), !Cluster %in% c("Normal", "Cortical admixture", "C2C")) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="PDL1", fill="Cluster")) + geom_boxplot(show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)[-(1:2)]] %>% unname)
print(p + xlab("") + ylab("CD274 RPPA (TCGA)") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
p = ggplot(bulk.meta %>% mutate(PDL1=bulk.mat2["CD274",] %>% as.numeric) %>% filter(!is.na(PDL1)) %>% filter(!Cluster  %in% c("Normal", "Cortical admixture", "C2C")) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="PDL1", fill="Cluster")) + geom_boxplot(show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)[-(1:2)]] %>% unname)
print(p + xlab("") + ylab("CD274 expression") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
dev.off()

pdf(file.path(fig3, "D_tcga_pdl1_rppa_and_mrna_violin.pdf"))
p = ggplot(bulk.meta %>% mutate(PDL1=tcga.rppa[f,"PDL1"]) %>% filter(!is.na(PDL1)) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="PDL1", fill="Cluster")) +  geom_violin(show.legend = F) + geom_point(alpha=0.5, position=position_jitter(width = 0.1), show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)[-1]] %>% unname)
print(p + xlab("") + ylab("CD274 RPPA (TCGA)") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
p = ggplot(bulk.meta %>% mutate(PDL1=bulk.mat2["CD274",] %>% as.numeric) %>% filter(!is.na(PDL1)) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="PDL1", fill="Cluster")) + geom_violin(show.legend = F) + geom_point(alpha=0.5, position=position_jitter(width = 0.1), show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)] %>% unname)
print(p + xlab("") + ylab("CD274 expression") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
dev.off()
pdf(file.path(fig3, "D_tcga_pdl1_rppa_and_mrna_violin_no_normals.pdf"))
p = ggplot(bulk.meta %>% mutate(PDL1=tcga.rppa[f,"PDL1"]) %>% filter(!is.na(PDL1), !Cluster %in% c("Normal", "Cortical admixture", "C2C")) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="PDL1", fill="Cluster")) +  geom_violin(show.legend = F) + geom_point(alpha=0.5, position=position_jitter(width = 0.1), show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)[-(1:2)]] %>% unname)
print(p + xlab("") + ylab("CD274 RPPA (TCGA)") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
p = ggplot(bulk.meta %>% mutate(PDL1=bulk.mat2["CD274",] %>% as.numeric) %>% filter(!is.na(PDL1)) %>% filter(!Cluster  %in% c("Normal", "Cortical admixture", "C2C")) %>% mutate(Cluster = add.counts(Cluster)), aes_string(x="Cluster", y="PDL1", fill="Cluster")) + geom_violin(show.legend = F) + geom_point(alpha=0.5, position=position_jitter(width = 0.1), show.legend = F) + scale_fill_manual(values=cluster.cols[levels(bulk.meta$Cluster)[-(1:2)]] %>% unname)
print(p + xlab("") + ylab("CD274 expression") + theme(text = element_text(size=16), axis.text.x = element_text(angle=-45, vjust=1, hjust=0), plot.margin = margin(l=1, t=4,r=3,unit="cm")))
dev.off()

# Brady Bunch

x = bulk.meta
f=match(bulk.meta$Sample, substr(tcga.rppa$Sample_ID,1,15) %>% make.names)
x$PDL1_rppa = tcga.rppa[f,"PDL1"]
x$PDL1_rna = as.numeric(bulk.mat2["CD274",])
x[,names(pseudobulk.normal.signatures)] = t(bulk.gsva)[,names(pseudobulk.normal.signatures)]

cols = setdiff(levels(sn.PPGLs$Cell_Type), "Tumour") %>% c("PDL1_rppa", "PDL1_rna", "Leukocyte_score")
labels = setdiff(levels(sn.PPGLs$Cell_Type), "Tumour") %>% c("CD274 RPPA (TCGA)", "CD274 expression", "Leukocyte score")

#p = lapply(1:length(cols), function(i){
#  p = ggplot(x, aes_string(x=make.names(cols[i]), y="Cluster", fill="Cluster")) + geom_hline(yintercept=1:length(unique(x[,"Cluster"]))) + geom_density_ridges(show.legend = F) + cluster.scale2 + scale_y_discrete(breaks=rev(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster)))) + theme(axis.ticks.y = element_blank()) + ylab("") + xlab(labels[i])
#  if(grepl("GSVA", labels[i])){
#    p = p + xlim(-1,1)
#  } 
#  return(p)
#})
#p[["ncol"]] = 3
#library(gridExtra)

labels=factor(labels, levels=labels)
p = ggplot()
for(i in 1:12){
  #p = p + geom_hline(data=data.frame(x,facet=labels[i]), yintercept = (1:length(unique(x$Cluster))))
  p = p + geom_density_ridges(data=data.frame(x, facet=labels[i]), mapping=aes_(y=as.integer(x$Cluster),x=as.name(ifelse(i <= 9, make.names(cols[i]), cols[i])), fill=~Cluster), show.legend=i==1, size=0.25, linetype='solid')
  if(i <= 9){
    p = p + geom_point(data=data.frame(facet=labels[i], x=c(-1,1), y=as.integer(1)), mapping=aes(x=x,y=y), alpha=0) # hack in [-1,1] range
  }
}
p = p + facet_wrap(~facet, ncol=3, scales="free_x") + cluster.scale2 + theme(panel.grid.major.y = element_line()) + scale_y_continuous(breaks=1:length(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster))), trans="reverse") + theme(axis.ticks.y = element_blank()) + ylab("") + xlab("") + theme(strip.background = element_blank()) + t + guides(fill = guide_legend(override.aes = list(shape = 15, size=6, col=NA))) #+ scale_x_continuous(limits=c(-1,1))) 

pdf(file.path(fig3, "Brady_bunch_ridge_plots.pdf"), width=10)
print(p)
dev.off()

labels=factor(labels, levels=labels)

pdf(file.path(fig3, "Brady_bunch_ridge_plots_just_gsva.pdf"), width=7)
for(idx in list(1:9, (1:5), 6:9)){
p = ggplot()
for(i in idx){
  p = p + geom_density_ridges(data=data.frame(x, facet=labels[i]), mapping=aes_(y=as.integer(x$Cluster),x=as.name(ifelse(i <= 9, make.names(cols[i]), cols[i])), fill=~Cluster), show.legend=i==idx[1], size=0.25, linetype='solid')
  if(i <= 9){
    p = p + geom_point(data=data.frame(facet=labels[i], x=c(-1,1), y=as.integer(1)), mapping=aes(x=x,y=y), alpha=0) # hack in [-1,1] range
  }
}
p = p + facet_wrap(~facet, ncol=round(sqrt(length(idx))), scales="fixed") + cluster.scale2 + theme(panel.grid.major.y = element_line()) + scale_y_continuous(breaks=1:length(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster))), trans="reverse") + theme(axis.ticks.y = element_blank()) + ylab("") + xlab("") + theme(strip.background = element_blank()) + t + guides(fill = guide_legend(override.aes = list(shape = 15, size=6, col=NA)))  #+ scale_x_continuous(limits=c(-1,1))) 
p = p + coord_equal(ratio=1/6)
print(p)
}
dev.off()

#p = lapply(1:length(cols), function(i){
#  p=ggplot(x, aes_string(y=make.names(cols[i]), x="Cluster", fill="Cluster")) + geom_boxplot(show.legend = F) + cluster.scale2 + scale_x_discrete(breaks=(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster)))) + theme(axis.ticks.x = element_blank()) + xlab("") + ylab(labels[i]) + theme(text = element_text(size=8))
#  if(grepl("GSVA", labels[i])){
#    p = p + ylim(-1,1)
#  } 
#  return(p)
#})
#p[["ncol"]] = 3
p = ggplot()
for(i in 1:12){
  p = p + geom_boxplot(data=data.frame(x, facet=labels[i]), mapping=aes_(x=as.integer(x$Cluster),y=as.name(ifelse(i <= 9, make.names(cols[i]), cols[i])), fill=~Cluster), show.legend=i==1)
  if(i <= 9){
    p = p + geom_point(data=data.frame(facet=labels[i], y=c(-1,1), x=as.integer(1)), mapping=aes(x=x,y=y), alpha=0) # hack in [-1,1] range
  }
}
p = p + facet_wrap(~facet, ncol=3, scales="free_y") + cluster.scale2 + theme(panel.grid.major.y = element_blank()) + scale_x_continuous(breaks=1:length(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster)))) + theme(axis.ticks.x = element_blank()) + ylab("") + xlab("") + theme(strip.background = element_blank()) + t + guides(fill = guide_legend(override.aes = list(shape = 15, size=6, col=NA))) #+ scale_x_continuous(limits=c(-1,1))) 


pdf(file.path(fig3, "Brady_bunch_boxplots.pdf"), width=10)
#do.call(grid.arrange, p)
p
dev.off()

#p = lapply(1:length(cols), function(i){
#  p=ggplot(x, aes_string(y=make.names(cols[i]), x="Cluster", fill="Cluster")) + geom_boxplot(show.legend = F) + cluster.scale2 + scale_x_discrete(breaks=rev(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster)))) + theme(axis.ticks.y = element_blank()) + xlab("") + ylab(labels[i]) + coord_flip()
#  if(grepl("GSVA", labels[i])){
#    p = p + ylim(-1,1)
#  } 
#  return(p)
#})
#p[["ncol"]] = 3

pdf(file.path(fig3, "Brady_bunch_boxplots_horizontal.pdf"), width=10)
#do.call(grid.arrange, p)
p + facet_wrap(~facet, ncol=3, scales="free_x") + coord_flip() + theme(axis.ticks.y = element_blank(), axis.ticks.x = element_line()) + scale_x_continuous(breaks=1:length(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster))), trans='reverse')
dev.off()

#genes = c("PNMT", "DDC", "DBH", "TH")
#p = lapply(genes, function(gene){
#  x[,gene] = as.numeric(bulk.mat2[gene,])
#  p=ggplot(x, aes_string(x=gene, y="Cluster", fill="Cluster"))  + geom_hline(yintercept=1:length(unique(x[,"Cluster"]))) + geom_density_ridges(show.legend = F) + cluster.scale2 + scale_y_discrete(breaks=rev(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster)))) + theme(axis.ticks.y = element_blank()) + ylab("") + xlab(gene)
#  if(grepl("GSVA", labels[i])){
#    p = p + ylim(-1,1)
#  } 
#  return(p)
#})
#p[["ncol"]] = 2
#pdf(file.path(fig3, "Brady_bunch_chromaffin_marker_ridge_plots.pdf"), width=7)
#do.call(grid.arrange, p)
#dev.off()

genes = c("PNMT", "DDC", "DBH", "TH")
x[,genes] = t(bulk.mat2[genes,])
x.range = range(bulk.mat2[genes,])
labels = factor(genes, levels=genes)

p = ggplot()
for(i in 1:length(genes)){
  gene=genes[i]
  p = p + geom_density_ridges(show.legend=i==1, data=data.frame(x, facet=labels[i], gene=x[,gene]), mapping=aes(y=as.integer(Cluster),fill=Cluster, x=gene), size=0.25)
}
pdf(file.path(fig3, "Brady_bunch_chromaffin_marker_ridge_plots.pdf"))
p + facet_wrap(~facet, ncol=2, scales="free_x", strip.position = "top") + cluster.scale2 + theme(strip.background = element_blank(), strip.text = element_text(face = 'italic'), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab("") + xlab("") + 
  theme(legend.background = element_blank(),
            legend.spacing = unit(0.1, "mm"),
            legend.key.size = unit(2, "mm")) + 
  guides(fill = guide_legend(override.aes = list(shape = 15, size=6, col=NA))) +
  theme(panel.grid.major.y = element_line()) + scale_y_continuous(breaks=1:length(levels(x$Cluster)), labels=(levels(x$Cluster)), trans='reverse')
dev.off()

### No normals#####

no.normals=!bulk.meta$Cluster %in% c("Normal", "C2C")
x = bulk.meta
f=match(bulk.meta$Sample, substr(tcga.rppa$Sample_ID,1,15) %>% make.names)
x$PDL1_rppa = tcga.rppa[f,"PDL1"]
x$PDL1_rna = as.numeric(bulk.mat2["CD274",])
x[,names(pseudobulk.normal.signatures)] = t(bulk.gsva)[,names(pseudobulk.normal.signatures)]
x = x %>% filter(no.normals)
x$Cluster = factor(x$Cluster %>% as.character, levels=intersect(levels(x$Cluster), unique(x$Cluster)))
cols = setdiff(levels(sn.PPGLs$Cell_Type), "Tumour") %>% c("PDL1_rppa", "PDL1_rna", "Leukocyte_score")
labels = setdiff(levels(sn.PPGLs$Cell_Type), "Tumour") %>% c("CD274 RPPA (TCGA)", "CD274 expression", "Leukocyte score")

labels=factor(labels, levels=labels)
p = ggplot()
for(i in 1:12){
  #p = p + geom_hline(data=data.frame(x,facet=labels[i]), yintercept = (1:length(unique(x$Cluster))))
  p = p + geom_density_ridges(data=data.frame(x, facet=labels[i]), mapping=aes_(y=as.integer(x$Cluster),x=as.name(ifelse(i <= 9, make.names(cols[i]), cols[i])), fill=~Cluster), show.legend=i==1, size=0.25, linetype='solid')
  if(i <= 9){
    p = p + geom_point(data=data.frame(facet=labels[i], x=c(-1,1), y=as.integer(1)), mapping=aes(x=x,y=y), alpha=0) # hack in [-1,1] range
  }
}
p = p + facet_wrap(~facet, ncol=3, scales="free_x") + cluster.scale2 + theme(panel.grid.major.y = element_line()) + scale_y_continuous(breaks=1:length(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster))), trans="reverse") + theme(axis.ticks.y = element_blank()) + ylab("") + xlab("") + theme(strip.background = element_blank()) + t + guides(fill = guide_legend(override.aes = list(shape = 15, size=6, col=NA))) #+ scale_x_continuous(limits=c(-1,1))) 

pdf(file.path(fig3, "Brady_bunch_ridge_plots_no_normals.pdf"), width=10)
print(p)
dev.off()

p = ggplot()
for(i in 1:12){
  p = p + geom_boxplot(data=data.frame(x, facet=labels[i]), mapping=aes_(x=as.integer(x$Cluster),y=as.name(ifelse(i <= 9, make.names(cols[i]), cols[i])), fill=~Cluster), show.legend=i==1)
  if(i <= 9){
    p = p + geom_point(data=data.frame(facet=labels[i], y=c(-1,1), x=as.integer(1)), mapping=aes(x=x,y=y), alpha=0) # hack in [-1,1] range
  }
}
p = p + facet_wrap(~facet, ncol=3, scales="free_y") + cluster.scale2 + theme(panel.grid.major.y = element_blank()) + scale_x_continuous(breaks=1:length(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster)))) + theme(axis.ticks.x = element_blank()) + ylab("") + xlab("") + theme(strip.background = element_blank()) + t + guides(fill = guide_legend(override.aes = list(shape = 15, size=6, col=NA))) #+ scale_x_continuous(limits=c(-1,1))) 


pdf(file.path(fig3, "Brady_bunch_boxplots_no_normals.pdf"), width=10)
#do.call(grid.arrange, p)
p
dev.off()

pdf(file.path(fig3, "Brady_bunch_boxplots_horizontal_no_normals.pdf"), width=10)
#do.call(grid.arrange, p)
p + facet_wrap(~facet, ncol=3, scales="free_x") + coord_flip() + theme(axis.ticks.y = element_blank(), axis.ticks.x = element_line()) + scale_x_continuous(breaks=1:length(levels(bulk.meta$Cluster)), labels=rep("", length(levels(bulk.meta$Cluster))), trans='reverse')
dev.off()

genes = c("PNMT", "DDC", "DBH", "TH")
x[,genes] = t(bulk.mat2[genes,no.normals])
x.range = range(bulk.mat2[genes,no.normals])
labels = factor(genes, levels=genes)

p = ggplot()
for(i in 1:length(genes)){
  gene=genes[i]
  p = p + geom_density_ridges(show.legend=i==1, data=data.frame(x, facet=labels[i], gene=x[,gene]), mapping=aes(y=as.integer(Cluster),fill=Cluster, x=gene), size=0.25)
}
pdf(file.path(fig3, "Brady_bunch_chromaffin_marker_ridge_plots_no_normals.pdf"))
p + facet_wrap(~facet, ncol=2, scales="free_x", strip.position = "top") + cluster.scale2 + theme(strip.background = element_blank(), strip.text = element_text(face = 'italic'), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + ylab("") + xlab("") + 
  theme(legend.background = element_blank(),
        legend.spacing = unit(0.1, "mm"),
        legend.key.size = unit(2, "mm")) + 
  guides(fill = guide_legend(override.aes = list(shape = 15, size=6, col=NA))) +
  theme(panel.grid.major.y = element_line()) + scale_y_continuous(breaks=1:length(levels(x$Cluster)), labels=(levels(x$Cluster)), trans='reverse')
dev.off()





# D T/NK cell umap

x = cbind(sn.PPGLs@meta.data, Embeddings(sn.PPGLs, reduction="T.cells"))
x=x[!is.na(x[,"T_1"]),]
t.umap.ratio = square.ratio(x[,c("T_1", "T_2")])
genes = "IFNG, CD2, GNLY, NKG7, GZMB, CD8A, CD4, FOXP3, PDCD1, CTLA4" %>% strsplit(split = ', ', fixed=T) %>% unlist
x[,genes] = sn.PPGLs@assays$SCT@data[genes,rownames(x)] %>% t
x$scMatch = ifelse(grepl("T cell|killer", x$scMatch), x$scMatch, NA)
labels = sapply(c("Cluster", "scMatch", genes), function(x){if(!x %in% genes){x}else{bquote(italic(.(x)))}})
labels = factor(labels, levels=labels)
pt.size=0.2
p = ggplot()
p = p + geom_point(data=data.frame(x, facet=factor("Cluster", levels=levels(labels))), mapping=aes(x=T_1,y=T_2, col=T.cells_Cluster), size=pt.size) + t.umap.ratio + scale_color_discrete(name='Cluster') + xlab("UMAP_1") + ylab("UMAP_2")
p = p  + guides(colour_new = guide_legend(override.aes = list(shape = 15, size=6))) + new_scale_colour()
p = p + geom_point(data=data.frame(x, facet=factor("scMatch", levels=levels(labels))) %>% filter(!is.na(scMatch)), mapping=aes(x=T_1,y=T_2, col=scMatch), size=pt.size) + scale_color_manual(name="scMatch", values=brewer.pal(n=4,name="Set1"))
p = p  + guides(colour_new_new = guide_legend(override.aes = list(shape = 15, size=6))) + new_scale_colour()

for(i in 1:length(genes)){
  p = p + stat_summary_hex(data=data.frame(x, facet=labels[i+2]), mapping=aes_(x=~T_1, y=~T_2, z=as.name(genes[i])), bins=30, fun = "mean")
}
p = p + scale_fill_gradientn(name="Mean log expression", colours=c("lightgrey", "red"))
p = p + facet_wrap(~facet, ncol=4, labeller = label_parsed)
p = p + theme(strip.background = element_blank())
pdf(file.path(fig3, "T_cell_marker_grid.pdf"))
print(p)
dev.off()

# Cell subtype dot plot



#################################################################
########################## Figure 4 #############################
#################################################################

fig4 = "Raw figures/Figure 4/Components"
dir.create(fig4, recursive=T, showWarnings = F)

# D dotplot
#genes = read.delim("Raw tables/Bulk/PPGLGenes.csv", sep=',', header=T, as.is=T)
#genes = read.delim("Raw tables/Bulk/PPGLGenes.tsv", sep='\t', header=T, as.is=T)
#genes = genes %>% filter(Select.for.dot.plot != "")
#genes = genes[,c("gene", "Select.for.dot.plot")]
#genes[,1] = gsub("^OPRM$", "OPRM1", genes[,1])
#genes=genes[!duplicated(genes),]
#genes = read.delim("Raw tables/Bulk/Curated genes for dotplot.txt", sep='\t', header=T, as.is=T)
genes = dotplot.genes[dotplot.genes$Figure=="4A",c("gene", "Pathway.x", "Order.x")]
genes[,1] = gsub("^OPRM$", "OPRM1", genes[,1])
genes[,2] = tools::toTitleCase(gsub("ransporter$", "ransporters", gsub("gudiance", "guidance", gsub("igaling", "ignaling", genes[,2]))))
genes[,2] = factor(genes[,2], levels=unique(genes[order(genes$Order.x),2]))

x = do.call(cbind, lapply(unique(sn.PPGLs$Sample), function(s){rowSums(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Sample==s & ((sn.PPGLs$Cell_Type %in% c("Chromaffin cells") & grepl("NAM", sn.PPGLs$orig.ident))|(sn.PPGLs$Cell_Type=="Tumour"))])})) %>%
  cbind(do.call(cbind, lapply(levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")), function(cell_type){
    rowSums(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Cell_Type==cell_type])
  }))) %>% as.matrix %>% LogNormalize(scale.factor=1e6) # logCPM
colnames(x) = c(unique(sn.PPGLs$Sample), levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")))
y = do.call(cbind, lapply(unique(sn.PPGLs$Sample), function(s){rowMeans(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Sample==s & ((sn.PPGLs$Cell_Type %in% c("Chromaffin cells") & grepl("NAM", sn.PPGLs$orig.ident))|(sn.PPGLs$Cell_Type=="Tumour"))] > 0)})) %>%
  cbind(do.call(cbind, lapply(levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")), function(cell_type){
    rowMeans(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Cell_Type==cell_type] > 0)
  }))) %>% as.matrix# logCPM
colnames(y) = c(unique(sn.PPGLs$Sample), levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")))
z = data.frame(Sample=colnames(y), stringsAsFactors = F)
z$Cell_Type = ifelse(z$Sample %in% sn.PPGLs$Cell_Type, z$Sample, ifelse(z$Sample %in% sn.PPGLs$Sample[which(sn.PPGLs$Genotype=="Normal")], "Chromaffin cells", "Tumour"))
f = c(which(z$Cell_Type=="Tumour")[order(sn.PPGLs$Genotype2[match(z$Sample[which(z$Cell_Type=="Tumour")], sn.PPGLs$Sample)])],
      which(z$Cell_Type=="Chromaffin cells"),
      which(z$Cell_Type %in% setdiff(levels(sn.PPGLs$Cell_Type), c("Tumour", "Chromaffin cells")))[order(sn.PPGLs$Cell_Type[match(z$Sample[which(z$Cell_Type %in% setdiff(levels(sn.PPGLs$Cell_Type), c("Tumour", "Chromaffin cells")))], sn.PPGLs$Cell_Type)])])
z=z[f,]
x=x[,f]
y=y[,f]
z$Genotype = sn.PPGLs$Genotype[match(z$Sample, sn.PPGLs$Sample)]
z$Cluster = sn.PPGLs$Subgroup4[match(z$Sample, sn.PPGLs$Sample)]



x=x[match(genes[,1], rownames(x)),]
y=y[match(genes[,1], rownames(y)),]
z[,"Cell Type"] = z$`Cell_Type`


pdf(file.path(fig4, "A_dotplots.pdf"), width=16, height=25)
hm.annot = HeatmapAnnotation(df = z[,c("Genotype", "Cluster", "Cell Type")], which="column", col = list("Genotype" = genotype.cols, "Cluster" = cluster.cols, "Cell Type" = cell.cols), show_legend = TRUE)
hm1 = HeatmapDotPlot(colour = x[,,drop=F], size = y[,,drop=F],
                       scale = TRUE, cell.size = 0.5,
                       cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
                       col=colorRamp2(c(-3,0,3),c("blue","grey90","red")),
                       column_split = (z$Cluster %>% (function(x){y=as.character(x); y[is.na(y)] = "Other cells"; factor(y, levels=c(levels(x), "Other cells"))})),
                       name="Scaled expression", row_names_gp=gpar(fontface="italic"),
                       row_split = genes[,2],
                       row_title_rot = 0)
draw(hm1)
#for(pathway in unique(genes[,2])){
#  print(pathway)
#  f2 = genes[,2]==pathway
#  hm1 = HeatmapDotPlot(colour = x[f2,,drop=F], size = y[f2,,drop=F],
#                       scale = TRUE, cell.size = 0.5,
#                       cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
#                       col=colorRamp2(c(-3,0,3),c("blue","white","red")),
#                       column_split = (z$Cluster %>% (function(x){y=as.character(x); y[is.na(y)] = "Other cells"; factor(y, levels=c(levels(x), "Other cells"))})),
#                       name="Scaled expression", row_names_gp=gpar(fontface="italic"),
#                       row_split = factor(genes[,2], levels=sort(unique(genes[,2])))[f2],
#                       row_title_rot = 0)
#  draw(hm1)
#}
dev.off()

pdf(file.path(fig4, "A_dotplots_no_normals.pdf"), width=16, height=25)
f=which(!is.na(z$Genotype))
hm.annot = HeatmapAnnotation(df = z[f,c("Genotype", "Cluster")], which="column", col = list("Genotype" = genotype.cols, "Cluster" = cluster.cols), show_legend = TRUE)
hm1 = HeatmapDotPlot(colour = x[,f,drop=F], size = y[,f,drop=F],
                     scale = TRUE, cell.size = 0.5,
                     cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
                     col=colorRamp2(c(-3,0,3),c("blue","grey90","red")),
                     column_split = (z$Cluster[f]),
                     name="Scaled expression", row_names_gp=gpar(fontface="italic"),
                     row_split = genes[,2],
                     row_title_rot = 0)
draw(hm1)


#f=which(!is.na(z$Genotype))
#hm.annot = HeatmapAnnotation(df = z[f,c("Genotype", "Cluster")], which="column", col = list("Genotype" = genotype.cols, "Cluster" = cluster.cols), show_legend = TRUE)
#for(pathway in unique(genes[,2])){
#  print(pathway)
#  f2 = genes[,2]==pathway
#  hm1 = HeatmapDotPlot(colour = x[f2,f,drop=F], size = y[f2,f,drop=F],
#                       scale = TRUE, cell.size = 0.5,
#                       cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
#                       col=colorRamp2(c(-3,0,3),c("blue","white","red")),
#                       column_split = z$Cluster[f],
#                       name="Scaled expression", row_names_gp=gpar(fontface="italic"),
#                       row_split = factor(genes[,2], levels=sort(unique(genes[,2])))[f2],
#                       row_title_rot = 0)
#  draw(hm1)
#}
dev.off()


#f2 = bulk.ccp[[9]]$consensusTree$order
#f3 = f2[order(factor(as.character(bulk.meta$Cluster[f2]), levels=c(setdiff(levels(bulk.meta$Cluster), c("Cortical admixture", "Normal", "C2C")), "C2C", "Normal")))]
#hm.bulk.annot = HeatmapAnnotation(df = bulk.meta[f3,c("Genotype", "Cluster")], col =list("Genotype" = genotype.cols, Cluster = cluster.cols), show_annotation_name = F)
#hm2 = Heatmap(bulk.mat2[match(genes[,1], rownames(bulk.mat2)),f3] %>% t %>% scale %>% t, bottom_annotation=hm.bulk.annot, col = colorRamp2(c(-3,0,3),c("blue","white","red")), show_heatmap_legend = FALSE, cluster_columns = F, cluster_rows=FALSE, show_column_names = FALSE, show_row_names = F, row_names_gp = gpar(fontface="italic"))






# Cell type corrplot
library(reshape)
f=which(sn.PPGLs$Genotype=="Normal")
x = table(sn.PPGLs$Sample[-f], gsub("Tumour", "Chromaffin cells", sn.PPGLs$Cell_Type[-f]))
x = x[!rownames(x) %in% c("Chromaffin cells", "Adrenocortical cells"),!colnames(x) %in% c("Chromaffin cells", "Adrenocortical cells")]
y1 = Hmisc::rcorr(x)
y2 = t(bulk.gsva[make.names(colnames(y1$r)),!bulk.meta$Cluster %in% c("Normal", "C2C")])
colnames(y2) = colnames(x)[match(colnames(y2), make.names(colnames(x)))]
y2 = Hmisc::rcorr(y2)
z = rbind(melt(y1$r) %>% mutate(Type="SN"), melt(y2$r) %>% mutate(Type="Bulk"))
colnames(z)[1:3] = c("x", "y", "cor")
z$p.value = c(melt(y1$P)[,3],melt(y2$P)[,3])
z$p.label = sapply(z$p.value, function(p){ifelse(p <= 0.05, paste(rep("*", min(3, as.integer(abs(log10(max(p,0.0001)))))), collapse=""), "")})
z[,1] = factor(z[,1], levels=intersect(levels(sn.PPGLs$Cell_Type), unique(z[,1])))
z[,2] = factor(z[,2], levels=intersect(levels(sn.PPGLs$Cell_Type), unique(z[,1])))
z = z %>% filter((Type=="Bulk" & as.integer(x) > as.integer(y)) | (Type=="SN" & as.integer(x) < as.integer(y)))
p = ggplot(z, aes(x=x,y=y, fill=cor)) + geom_tile(width=0.9,height=0.9) + scale_fill_gradientn(colours=c("blue","white","red"), na.value = "white", name=expression(italic("Rho")), limits=c(-1,1)) + geom_text(mapping = aes(x=x,y=y,label=round(cor, digits=2))) + geom_text(mapping = aes(x=x,y=as.integer(y)+0.2,label=p.label)) + geom_abline(slope=1,intercept=0, col="grey50")
pdf(file.path(fig4, "cell_type_correlation.pdf"))
p + theme_minimal() + theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 90, hjust=1)) + xlab("Bulk") + ylab("Single-nuclei") + coord_equal()
dev.off()

# NATMI counts heatmap
natmi.cols = c(col1 = "dodgerblue4", col2 = 'peachpuff', col3 = 'deeppink4')

# summed specificity?

natmi.results$Sending.cluster = gsub("Macrophages.Monocytes", "Myeloid cells", gsub("T cells", "T/NK cells", natmi.results$Sending.cluster))
natmi.results$Target.cluster =gsub("Macrophages.Monocytes", "Myeloid cells", gsub("T cells", "T/NK cells", natmi.results$Target.cluster))
natmi.results$Sending.cluster = factor(natmi.results$Sending.cluster, levels = levels(sn.PPGLs$Cell_Type))
natmi.results$Target.cluster = factor(natmi.results$Target.cluster, levels = levels(sn.PPGLs$Cell_Type))
natmi.results$Subgroup = gsub("EPAS1", "PH-NOS", gsub("MAML$", "MAML3", natmi.results$Subgroup))
#x = do.call(cbind, lapply(levels(sn.PPGLs$Subgroup4) %>% setdiff("C2C"), function(cluster){
#  print(cluster)
#  cluster = sn.PPGLs$Subgroup3[match(cluster, sn.PPGLs$Subgroup4)]
x = do.call(cbind, lapply(c("Tumour", "Normal"), function(type){
  if(type=="Tumour"){cluster = setdiff(unique(natmi.results$Subgroup), "Normal")} else {cluster="Normal"}
  x = natmi.results %>% filter(Subgroup%in%cluster) %>%
    filter(Ligand.detection.rate >= 0.2, Receptor.detection.rate >= 0.2, Ligand.expressing.cells/Ligand.detection.rate >= 30,
           Receptor.expressing.cells/Receptor.detection.rate >= 30) %>%
    filter(!((Sending.cluster %in% c("Chromaffin cells") | Target.cluster == "Chromaffin cells") & Genotype != "Normal")) %>%
    mutate(Sending.cluster = factor(gsub("Tumour", "Chromaffin cells", Sending.cluster), levels=levels(Sending.cluster)),
           Target.cluster = factor(gsub("Tumour", "Chromaffin cells", Target.cluster), levels=levels(Target.cluster)))
  
  y = x %>% group_by(Sending.cluster, Target.cluster, orig.ident) %>% summarise(x=sum(Edge.average.expression.derived.specificity)) %>%
    group_by(Sending.cluster, Target.cluster) %>% summarise(x=mean(x))
  m = sparseMatrix(i = as.integer(y$Sending.cluster), j=as.integer(y$Target.cluster),
                   x=y$x, dims=c(length(levels(y$Sending.cluster)), length(levels(y$Target.cluster))),
                   dimnames = list(levels(y$Sending.cluster), levels(y$Target.cluster))) %>% as.matrix()
  m = m[!grepl("Tumour", rownames(m)), !grepl("Tumour", colnames(m))]
  return(m)
}))

hm.annot = HeatmapAnnotation(df = data.frame(`Cell Type` = colnames(x), check.names=F), which="column", col = list("Cell Type" = cell.cols), show_annotation_name=F, show_legend=F)
hm.annot2 = HeatmapAnnotation(df = data.frame(`Cell Type` = rownames(x), check.names=F), which="row", col = list("Cell Type" = cell.cols), show_annotation_name = F, show_legend=F)
clusters = unlist(lapply(levels(sn.PPGLs$Subgroup4) %>% setdiff("C2C"), function(cluster) rep(cluster, length(unique(natmi.results$Sending.cluster))-1)))
f=!(sn.PPGLs$Genotype != "Normal" & sn.PPGLs$Cell_Type=="Chromaffin cells")
tab = table(sn.PPGLs$orig.ident[f], gsub("Tumour", "Chromaffin cells", sn.PPGLs$Cell_Type[f]))
hm.annot3 = HeatmapAnnotation(`Cell Proportion` = anno_boxplot(x = lapply(1:ncol(x), function(i){
  cluster = clusters[i]
  cell_type = colnames(x)[i]
  if(cluster != "Normal" & cell_type == "Chromaffin cell"){cell_type = "Tumour"}
  f=rownames(tab) %in% unique(sn.PPGLs$orig.ident[sn.PPGLs$Subgroup4 == cluster])
  return(tab[f,cell_type]/rowSums(tab)[f])
}), gp=gpar(fill=cell.cols[colnames(x)])), which="column")

clusters = c(rep("Tumour", 9), rep("Normal", 9))
rownames(x) = gsub("^Chromaffin", "Tumour/Chromaffin", rownames(x))
colnames(x) = sapply(colnames(x), function(x){
  if(grepl("Myeloid", x)){return("My")}
  substr(x,1,1)
})
hm = Heatmap(x, cluster_rows = F, cluster_columns=F, #col=colorRamp2(breaks = c(0, max(x)/2, max(x)), natmi.cols),
             name="Average summed specificity weight", height=unit(0.5*nrow(x), "cm"), width=unit(0.5*ncol(x),"cm"),
             column_split = clusters,
             show_column_names=T, bottom_annotation = hm.annot, right_annotation = hm.annot2,
             row_title="Ligand expressing", row_title_side="left", column_title_side="top")
pdf(file.path(fig4, "NATMI_heatmap_tumour_vs_normal.pdf"))
draw(hm, heatmap_legend_side = "top")
decorate_heatmap_body("Average summed specificity weight", {
  grid.text("Receptor expressing", just=c("center", "bottom"), y=unit(0, "npc") - unit(cell.size*3.5, "cm"),x=unit(0.5 + 0.5, "npc"), gp = gpar(fontsize = 14))
})
dev.off()

# add natmi annotation
genes = read.delim("Raw tables/Bulk/PPGLGenes.tsv", sep='\t', header=T, as.is=T)
genes$Natmi.Ligand = genes$gene %in% natmi.results$Ligand.symbol
genes$Natmi.Receptor = genes$gene %in% natmi.results$Receptor.symbol


# immune cell subtyping
GSE146771_10x = read.delim("data/GSE146771_CRC.Leukocyte.10x.TPM.txt", header=T, as.is=T, sep=' ', row.names=1, quote='"')
GSE146771_meta = read.delim("data/GSE146771_CRC.Leukocyte.10x.Metadata.txt", header=T, as.is=T, sep='\t', quote='"')
GSE146771_10x_centroids = do.call(cbind, lapply(unique(GSE146771_meta$Sub_Cluster) %>%(function(x){setNames(x,x)}), function(cluster){rowMeans(GSE146771_10x[,which(GSE146771_meta$Sub_Cluster==cluster)])}))

x = sn.PPGLs@assays$SCT@counts[,which(sn.PPGLs$Cell_Type %in% c("Mast cells", "Myeloid cells", "T/NK cells", "B cells"))]


myeloid.signatures = read.delim("Raw tables/Pseudobulk/myeloid_cell_signatures.txt", sep='\t', header=T, as.is=T, check.names=F, skip=2)
colnames(myeloid.signatures)[1:5] = c("Gene", "Cluster", "AUC", "F_statistic", "F_adjusted_pvalue")
myeloid.signatures=myeloid.signatures[!is.na(myeloid.signatures[,3]),!is.na(myeloid.signatures[1,])]
myeloid.signatures = myeloid.signatures[!duplicated(myeloid.signatures[,1]), -(2:5)]
rownames(myeloid.signatures) = myeloid.signatures[,1]
myeloid.signatures = myeloid.signatures[,-1]
myeloid.signatures = as.matrix(myeloid.signatures)
myeloid.signatures = myeloid.signatures[intersect(rownames(myeloid.signatures), rownames(sn.PPGLs@assays$RNA@data)),]
x = sn.PPGLs@assays$RNA@counts[match(intersect(rownames(myeloid.signatures), rownames(sn.PPGLs@assays$RNA@counts)), rownames(sn.PPGLs@assays$RNA@counts)),which(sn.PPGLs$Cell_Type=="Myeloid cells")]
x = as.matrix(x)
y = cor(x, myeloid.signatures, method="spearman")
z = data.frame(Barcode = rownames(y), Top_Cluster = sapply(1:nrow(y), function(i){colnames(y)[y[i,]==max(y[i,])]}), Correlation = apply(y,1,max))

x = cbind(sn.PPGLs@meta.data, Embeddings(sn.PPGLs,reduction='Macrophages.Monocytes'))
x = x[!is.na(x$Macrophages_1),]
x[,colnames(z)] = z[rownames(x),]
r = square.ratio(x[,c("Macrophages_1", "Macrophages_2")])
pdf(file.path(fig4, "Myeloid_subtypes.pdf"))
ggplot(x, aes(x=Macrophages_1, y=Macrophages_2, col=Top_Cluster)) + geom_point() + scale_colour_manual(values=c25, name="Cell Type") + t + t2 + r + xlab("UMAP_1") + ylab("UMAP_2")
FeaturePlot(sn.PPGLs, features=c("VCAN", "IFITM2", "FCGR3A", "CD14"), reduction="Macrophages.Monocytes")
dev.off()


##### Figure 5
fig5 = "Raw figures/Figure 5/Components"
if(!dir.exists(fig5)){dir.create(fig5, recursive = T)}
genes = dotplot.genes[dotplot.genes$Figure=="5A",c("gene", "Pathway.x", "Order.x")]
genes[,1] = gsub("^OPRM$", "OPRM1", genes[,1])
genes[,2] = tools::toTitleCase(gsub("ransporter$", "ransporters", gsub("gudiance", "guidance", gsub("igaling", "ignaling", genes[,2]))))
genes[,2] = factor(genes[,2], levels=unique(genes[order(genes$Order.x),2]))

x = do.call(cbind, lapply(unique(sn.PPGLs$Sample), function(s){rowSums(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Sample==s & ((sn.PPGLs$Cell_Type %in% c("Chromaffin cells") & grepl("NAM", sn.PPGLs$orig.ident))|(sn.PPGLs$Cell_Type=="Tumour"))])})) %>%
  cbind(do.call(cbind, lapply(levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")), function(cell_type){
    rowSums(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Cell_Type==cell_type])
  }))) %>% as.matrix %>% LogNormalize(scale.factor=1e6) # logCPM
colnames(x) = c(unique(sn.PPGLs$Sample), levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")))
y = do.call(cbind, lapply(unique(sn.PPGLs$Sample), function(s){rowMeans(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Sample==s & ((sn.PPGLs$Cell_Type %in% c("Chromaffin cells") & grepl("NAM", sn.PPGLs$orig.ident))|(sn.PPGLs$Cell_Type=="Tumour"))] > 0)})) %>%
  cbind(do.call(cbind, lapply(levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")), function(cell_type){
    rowMeans(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Cell_Type==cell_type] > 0)
  }))) %>% as.matrix# logCPM
colnames(y) = c(unique(sn.PPGLs$Sample), levels(sn.PPGLs$Cell_Type) %>% setdiff(c("Tumour", "Chromaffin cells")))
z = data.frame(Sample=colnames(y), stringsAsFactors = F)
z$Cell_Type = ifelse(z$Sample %in% sn.PPGLs$Cell_Type, z$Sample, ifelse(z$Sample %in% sn.PPGLs$Sample[which(sn.PPGLs$Genotype=="Normal")], "Chromaffin cells", "Tumour"))
f = c(which(z$Cell_Type=="Tumour")[order(sn.PPGLs$Genotype2[match(z$Sample[which(z$Cell_Type=="Tumour")], sn.PPGLs$Sample)])],
      which(z$Cell_Type=="Chromaffin cells"),
      which(z$Cell_Type %in% setdiff(levels(sn.PPGLs$Cell_Type), c("Tumour", "Chromaffin cells")))[order(sn.PPGLs$Cell_Type[match(z$Sample[which(z$Cell_Type %in% setdiff(levels(sn.PPGLs$Cell_Type), c("Tumour", "Chromaffin cells")))], sn.PPGLs$Cell_Type)])])
z=z[f,]
x=x[,f]
y=y[,f]
z$Genotype = sn.PPGLs$Genotype[match(z$Sample, sn.PPGLs$Sample)]
z$Cluster = sn.PPGLs$Subgroup4[match(z$Sample, sn.PPGLs$Sample)]



x=x[match(genes[,1], rownames(x)),]
y=y[match(genes[,1], rownames(y)),]
z[,"Cell Type"] = z$`Cell_Type`


pdf(file.path(fig4, "A_dotplots.pdf"), width=16, height=25)
hm.annot = HeatmapAnnotation(df = z[,c("Genotype", "Cluster", "Cell Type")], which="column", col = list("Genotype" = genotype.cols, "Cluster" = cluster.cols, "Cell Type" = cell.cols), show_legend = TRUE)
hm1 = HeatmapDotPlot(colour = x[,,drop=F], size = y[,,drop=F],
                     scale = TRUE, cell.size = 0.5,
                     cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
                     col=colorRamp2(c(-3,0,3),c("blue","grey90","red")),
                     column_split = (z$Cluster %>% (function(x){y=as.character(x); y[is.na(y)] = "Other cells"; factor(y, levels=c(levels(x), "Other cells"))})),
                     name="Scaled expression", row_names_gp=gpar(fontface="italic"),
                     row_split = genes[,2],
                     row_title_rot = 0)
draw(hm1)

dev.off()

pdf(file.path(fig4, "A_dotplots_no_normals.pdf"), width=16, height=25)
f=which(!is.na(z$Genotype))
hm.annot = HeatmapAnnotation(df = z[f,c("Genotype", "Cluster")], which="column", col = list("Genotype" = genotype.cols, "Cluster" = cluster.cols), show_legend = TRUE)
hm1 = HeatmapDotPlot(colour = x[,f,drop=F], size = y[,f,drop=F],
                     scale = TRUE, cell.size = 0.5,
                     cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
                     col=colorRamp2(c(-3,0,3),c("blue","grey90","red")),
                     column_split = (z$Cluster[f]),
                     name="Scaled expression", row_names_gp=gpar(fontface="italic"),
                     row_split = genes[,2],
                     row_title_rot = 0)
draw(hm1)

dev.off()

cell_proportion_dotplot = function(samples, cell_types, all_samples, groups=NULL, barplot=FALSE, annot.label = "Total Proportions", ...){
  if(is.character(all_samples[1])){all_samples = factor(all_samples, levels=unique(all_samples))}
  samples=factor(samples, levels=levels(all_samples))
  sample_counts = table(all_samples)
  cell_table = table(cell_types, samples)
  cell_proportions_total = sweep(cell_table,2,pmax(1,sample_counts), "/") %>% Matrix %>% as.matrix
  cell_proportions = sweep(cell_table,2,pmax(1,colSums(cell_table)), "/") %>%  Matrix %>% as.matrix
  
  if(barplot){
    args = list(which="column")
    args[[annot.label]] = anno_barplot(sapply(colnames(cell_proportions) %>% (function(x){setNames(x,x)}), function(sample){cell_proportions_total[,sample] %>% sum}))
    hm.annot = do.call(HeatmapAnnotation, args)
    if(!is.null(groups)){
      hm = HeatmapDotPlot(colour = cell_proportions, size=cell_proportions, top_annotation=hm.annot, scale=FALSE, col = colorRamp2(c(0,1), c("lightgrey","black")),  row_split = groups[match(rownames(cell_proportions), cell_types)], ...)
    } else {
      hm = HeatmapDotPlot(colour = cell_proportions, size=cell_proportions, top_annotation=hm.annot, scale=FALSE, col = colorRamp2(c(0,1), c("lightgrey","black")), ...)
    }
    
  } else {
    if(!is.null(groups)){
      hm = HeatmapDotPlot(colour = cell_proportions, size=cell_proportions, scale=FALSE, col = colorRamp2(c(0,1), c("lightgrey","black")), groups[match(rownames(cell_proportions), cell_types)], ...)
    } else {
      hm = HeatmapDotPlot(colour = cell_proportions, size=cell_proportions, scale=FALSE, col = colorRamp2(c(0,1), c("lightgrey","black")), ...)
    }
  }
  return(hm)
}

#hm.total = cell_proportion_dotplot(sn.PPGLs$Sample, sn.PPGLs$Cell_Type, sn.PPGLs$Sample[order(sn.PPGLs$Genotype)], barplot=FALSE, top_annotation=HeatmapAnnotation(`Total cells`=anno_barplot(table(sn.PPGLs$Sample)[unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)])] %>% as.numeric), which="column"), column_split = sn.PPGLs$Subgroup4[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)], bottom_annotation = HeatmapAnnotation(which="column", df = data.frame(Genotype = sn.PPGLs$Genotype[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)]), col = list("Genotype" = genotype.cols)), column_title_rot = 90, size.label="Proportion", colour.label=NULL)
#hm.total
#hm = cell_proportion_dotplot(a$Sample, a$Myeloid_Subtype, sn.PPGLs$Sample[order(sn.PPGLs$Genotype)], barplot=TRUE,column_split = sn.PPGLs$Subgroup4[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)], bottom_annotation = HeatmapAnnotation(which="column", df = data.frame(Genotype = sn.PPGLs$Genotype[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)]), col = list("Genotype" = genotype.cols), show_annotation_name = FALSE), column_title_rot = 90, annot.label="Myeloid cell proportions")
#hms

#hm.list = hm.total$object %v% hm$object
#hm.legend = list(Legend(title = "Proportion", at=c(0.25, 0.50, 0.75, 1.00), labels = c(0.25, 0.50, 0.75, 1.00) %>% as.character, size=unit.c(unit(sqrt(0.25)/2, "cm"), unit(sqrt(0.5)/2, "cm"), unit(sqrt(0.75)/2 , "cm"), unit(sqrt(1.0)/2, "cm")), type = "points", grid_height = unit(0.5,"cm"), grid_width=unit(0.5,"cm"), legend_height=unit(4, "cm"), background = NA, legend_gp = gpar(col=colorRamp2(c(0,1), c("lightgrey","black"))(c(0.25,0.5,0.75,1)))))

#pdf(file.path(fig3, "Cell_type_proportion_dotplot.pdf"), width=16, height=16)
#draw(hm.list, annotation_legend_list = hm.legend)
#dev.off()

#sn.PPGLs$GSE131907_supertype = GSE131907_annot$Cell_type.refined[match(sn.PPGLs$GSE131907_subtype, GSE131907_annot$Cell_subtype)]
#sn.PPGLs$GSE146771_supertype = GSE146771_meta$Global_Cluster[match(sn.PPGLs$GSE146771_subtype, GSE146771_meta$Sub_Cluster)]

#hm = cell_proportion_dotplot(sn.PPGLs$Sample, sn.PPGLs$GSE146771_subtype, sn.PPGLs$Sample[order(sn.PPGLs$Genotype)], groups=sn.PPGLs$GSE146771_supertype, barplot=TRUE,column_split = sn.PPGLs$Subgroup4[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)], bottom_annotation = HeatmapAnnotation(which="column", df = data.frame(Genotype = sn.PPGLs$Genotype[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)]), col = list("Genotype" = genotype.cols), show_annotation_name = FALSE), column_title_rot = 90, annot.label="Immune cell proportions", row_title_rot=0)
#hm2 = sn.PPGLs@meta.data[immune.cells,] %>% (function(sn.PPGLs){cell_proportion_dotplot(sn.PPGLs$Sample, sn.PPGLs$GSE131907_subtype, sn.PPGLs$Sample[order(sn.PPGLs$Genotype)], groups=sn.PPGLs$GSE131907_supertype, barplot=TRUE,column_split = sn.PPGLs$Subgroup4[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)], bottom_annotation = HeatmapAnnotation(which="column", df = data.frame(Genotype = sn.PPGLs$Genotype[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)]), col = list("Genotype" = genotype.cols), show_annotation_name = FALSE), column_title_rot = 90, annot.label="Immune cell proportions", row_title_rot=0)})
#hm3 = sn.PPGLs@meta.data[other.normals,] %>% (function(sn.PPGLs){cell_proportion_dotplot(sn.PPGLs$Sample, sn.PPGLs$GSE131907_subtype, sn.PPGLs$Sample[order(sn.PPGLs$Genotype)], groups=sn.PPGLs$GSE131907_supertype, barplot=TRUE,column_split = sn.PPGLs$Subgroup4[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)], bottom_annotation = HeatmapAnnotation(which="column", df = data.frame(Genotype = sn.PPGLs$Genotype[match(unique(sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]), sn.PPGLs$Sample)]), col = list("Genotype" = genotype.cols), show_annotation_name = FALSE), column_title_rot = 90, annot.label="Immune cell proportions", row_title_rot=0)})

#hm4 = Reduce(function(x,y){x %v% y}, lapply(sn.PPGLs$GSE131907_supertype %>% unique %>% na.omit, function(supertype){
#  all.samples = sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]
#  sn.PPGLs$Sample = factor(sn.PPGLs$Sample, levels=unique(all.samples))
#  col.split = sn.PPGLs$Subgroup4[match(levels(sn.PPGLs$Sample), sn.PPGLs$Sample)]
#  genotypes = sn.PPGLs$Genotype[match(levels(sn.PPGLs$Sample), sn.PPGLs$Sample)]
#  print(supertype)
#  (sn.PPGLs@meta.data[sn.PPGLs$GSE131907_supertype == supertype,] %>% (function(sn.PPGLs){cell_proportion_dotplot(sn.PPGLs$Sample, sn.PPGLs$GSE131907_subtype, all.samples, groups=NULL, barplot=TRUE,column_split = col.split, bottom_annotation = HeatmapAnnotation(which="column", df = data.frame(Genotype = genotypes), col = list("Genotype" = genotype.cols), show_annotation_name = FALSE), column_title_rot = 90, annot.label=paste(supertype, "proportions"), row_title_rot=0)}))$object
#}))
#draw(hm4)
#hm5 = Reduce(function(x,y){x %v% y}, lapply(sn.PPGLs$GSE146771_supertype %>% unique %>% na.omit, function(supertype){
#  all.samples = sn.PPGLs$Sample[order(sn.PPGLs$Genotype)]
#  sn.PPGLs$Sample = factor(sn.PPGLs$Sample, levels=unique(all.samples))
#  col.split = sn.PPGLs$Subgroup4[match(levels(sn.PPGLs$Sample), sn.PPGLs$Sample)]
#  genotypes = sn.PPGLs$Genotype[match(levels(sn.PPGLs$Sample), sn.PPGLs$Sample)]
#  print(supertype)
#  (sn.PPGLs@meta.data[sn.PPGLs$GSE146771_supertype == supertype,] %>% (function(sn.PPGLs){cell_proportion_dotplot(sn.PPGLs$Sample, sn.PPGLs$GSE146771_subtype, all.samples, groups=NULL, barplot=TRUE,column_split = col.split, bottom_annotation = HeatmapAnnotation(which="column", df = data.frame(Genotype = genotypes), col = list("Genotype" = genotype.cols), show_annotation_name = FALSE), column_title_rot = 90, annot.label=paste(supertype, "proportions"), row_title_rot=0)}))$object
#}))
#draw(hm5)
#pdf('~/a.pdf', height=30, width=10)
#draw(hm4)
#draw(hm5)
#dev.off()
sample.order = sn.PPGLs$Sample[order(sn.PPGLs$Genotype)] %>% unique
cell.order = levels(sn.PPGLs$Cell_Supertype) %>% setdiff("Mast cells")
max.size=0.5
hm = table(sn.PPGLs$Cell_Supertype, sn.PPGLs$Sample)[,sample.order] %>%
  (function(x){sweep(x,2,colSums(x), "/")[cell.order,]}) %>% (function(x){y = do.call(cbind, lapply(1:ncol(x), function(i){x[,i]})); dimnames(y) = dimnames(x); return(y)}) %>%
  (function(x){x=as.matrix(x); HeatmapDotPlot(x, x,
                                              col = function(x){colorRamp2(c(0,max.size),c("lightgrey", "black"))(x) %>% alpha(1)},
                                              size.is.col=TRUE, min.size=0, max.size=max.size, size.label="Cell proportion", scale=FALSE,
                                              top_annotation = HeatmapAnnotation(df = sn.PPGLs@meta.data[match(sample.order, sn.PPGLs$Sample),c("Genotype"), drop=F],
                                                                                    col = list("Genotype" = genotype.cols), which = "column",
                                                                                 show_legend = TRUE, show_annotation_name=FALSE),
                                              column_split = sn.PPGLs$Subgroup4[match(sample.order, sn.PPGLs$Sample)],
                                              column_title_rot = 90, row_title_rot=0,
                                              row_split = rep("Major cell types", length(cell.order)),
                                              show_row_names=TRUE, show_column_names=FALSE)})
hm$object = hm$object %v% Reduce(function(x,y){x %v% y}, lapply(c("Endothelial cells", "Fibroblasts", "Myeloid cells", "T/NK cells", "B cells"), function(cell_type){
  f = which(sn.PPGLs$Cell_Supertype==cell_type)
  cell_subtypes = sort(unique(sn.PPGLs$Cell_Type[f]))
  hm = table(sn.PPGLs$Cell_Type[f], sn.PPGLs$Sample[f])[,sample.order] %>%
    (function(x){x=sweep(x,2,colSums(x), "/")[cell_subtypes,]; x[is.na(x)]=0; x}) %>% (function(x){y = do.call(cbind, lapply(1:ncol(x), function(i){x[,i]})); dimnames(y) = dimnames(x); return(y)}) %>%
    (function(x){x=as.matrix(x); rownames(x) = gsub("pos ", "+ ", gsub("neg ", "- ", rownames(x))); HeatmapDotPlot(x, x,
                                                col = function(x){colorRamp2(c(0,max.size),c("lightgrey", "black"))(x) %>% alpha(1)},
                                                size.is.col=TRUE, min.size=0, max.size=max.size, size.label="Cell proportion", scale=FALSE,
                                                column_split = sn.PPGLs$Subgroup4[match(sample.order, sn.PPGLs$Sample)],
                                                column_title_rot = 90, row_title_rot=0,
                                                row_split = rep(cell_type, length(cell_subtypes)),
                                                show_column_names=cell_type == "B cells")})
  return(hm$object)
}))
#pdf(file.path(fig3, "Cell_subtype_proportions_dotplot.pdf"))
do.call(draw, hm)
#dev.off()

hm2 = Reduce(function(x,y){x %v% y}, lapply(c("Endothelial cells", "Fibroblasts", "Myeloid cells", "T/NK cells", "B cells"), function(cell_type){
  f = which(sn.PPGLs$Cell_Supertype==cell_type)
  cell_subtypes = sort(unique(sn.PPGLs$Cell_Type[f]))
  hm = table(sn.PPGLs$Cell_Type[f], sn.PPGLs$Sample[f])[,sample.order] %>%
    (function(x){x=sweep(x,2,colSums(x), "/")[cell_subtypes,]; x[is.na(x)]=0; x}) %>%
    (function(x){y = do.call(cbind, lapply(1:ncol(x), function(i){x[,i]})); dimnames(y) = dimnames(x); return(y)}) %>%
    (function(x){x=as.matrix(x); rownames(x) = gsub("pos ", "+ ", gsub("neg ", "- ", rownames(x)));
    print(cell_type)
    if(cell_type == 'B cells'){
      bottom_annot = columnAnnotation(Genotype = sn.PPGLs@meta.data[match(sample.order, sn.PPGLs$Sample),c("Genotype"),],
                                            col = list("Genotype" = genotype.cols),
                                            show_legend = TRUE, show_annotation_name=FALSE)
    } else {bottom_annot = NULL}
    HeatmapDotPlot(x, x,
                   col = function(x){colorRamp2(c(0,max.size),c("lightgrey", "black"))(x) %>% alpha(1)},
                   size.is.col=TRUE, min.size=0, max.size=max.size, size.label="Cell proportion", scale=FALSE,
                   column_split = sn.PPGLs$Subgroup4[match(sample.order, sn.PPGLs$Sample)],
                   column_title_rot = 90, row_title_rot=0,
                   row_split = rep(cell_type, length(cell_subtypes)),
                   show_column_names=cell_type == "B cells", bottom_annotation = bottom_annot,
                   top_annotation = HeatmapAnnotation(`Total Proportion` = anno_barplot(x = as.numeric(table(sn.PPGLs$Sample, sn.PPGLs$Cell_Supertype)[colnames(x),cell_type])/as.numeric(table(sn.PPGLs$Sample)[colnames(x)])), which='column'))})
  return(hm$object)
}))

pdf(file.path(fig3, "Cell_subtype_proportions_dotplot_with_total_proportion_barplots.pdf"), width=16, height=9)
draw(hm2)
dev.off()

# Supplementary
qc_metadata = read.delim("Raw tables/Pseudobulk/single_nuclei_QC_metadata.tsv", header=T, as.is=T, sep='\t')

supp1 = file.path("Raw figures", "Figure S1", "Components")
if(!dir.exists(supp1)){dir.create(supp1, recursive = T)}
pdf(file.path(supp1, "QC_cutoffs_boxplots.pdf"), width=10)

#p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=log2(nCount_RNA))) + geom_point(position=position_jitter(width=0.1), alpha=0.1, show.legend = F) + geom_violin(alpha=0.5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=log2(nCount_RNA))) + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p = p + geom_vline(aes(xintercept = factor(gsub("P018-PGL3", "P018-PGL1",Sample)) %>% as.integer() %>% (function(x) {x+0.5})), size=0.25)
p = p + xlab("") + ylab("Log2 UMI count")
print(p)

#p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=percent.mt))+ geom_point(position=position_jitter(width=0.1), alpha=0.1, show.legend = F) + geom_violin(alpha=0.5) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) 
p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=percent.mt))+ geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
p = p + geom_vline(aes(xintercept = factor(gsub("P018-PGL3", "P018-PGL1",Sample)) %>% as.integer() %>% (function(x) {x+0.5})), size=0.25)
p = p + xlab("") + ylab("Percent mitochondrial counts")
print(p)

#p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=nCount_RNA_MAD, fill=Immune_Cell, col=Immune_Cell)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_violin() + scale_fill_discrete(name='') + guides(col=FALSE)# + geom_point(position=position_jitterdodge(jitter.width=0.1), alpha=0.2)
p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=nCount_RNA_MAD, fill=Immune_Cell)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_boxplot() + scale_fill_discrete(name='') + guides(col=FALSE)# + geom_point(position=position_jitterdodge(jitter.width=0.1), alpha=0.2)
p = p + geom_hline(data = data.frame(Cutoff=c("Non-immune cell cutoff", "Immune cell cutoff") %>% (function(x){factor(x, levels=x)}), Value=c(-2.5, -4)), mapping=aes(yintercept=Value, linetype=Cutoff)) + scale_linetype_manual(values = c("Non-immune cell cutoff" = 5, "Immune cell cutoff" = 6), name='')
p = p + geom_vline(aes(xintercept = factor(gsub("P018-PGL3", "P018-PGL1",Sample)) %>% as.integer() %>% (function(x) {x+0.5})), size=0.25)
p = p + xlab("") + ylab("Log2 UMI count (MAD)")
print(p)

#p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=nFeature_RNA_MAD, fill=Immune_Cell, col=Immune_Cell)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_violin() + scale_fill_discrete(name='') + guides(col=FALSE)# + geom_point(position=position_jitterdodge(jitter.width=0.1), alpha=0.2)
p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=nFeature_RNA_MAD, fill=Immune_Cell)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_boxplot() + scale_fill_discrete(name='') + guides(col=FALSE)# + geom_point(position=position_jitterdodge(jitter.width=0.1), alpha=0.2)
p = p + geom_hline(data = data.frame(Cutoff=c("Non-immune cell cutoff", "Immune cell cutoff") %>% (function(x){factor(x, levels=x)}), Value=c(-2.5, -4)), mapping=aes(yintercept=Value, linetype=Cutoff)) + scale_linetype_manual(values = c("Non-immune cell cutoff" = 5, "Immune cell cutoff" = 6), name='')
p = p + geom_vline(aes(xintercept = factor(gsub("P018-PGL3", "P018-PGL1",Sample)) %>% as.integer() %>% (function(x) {x+0.5})), size=0.25)
p = p + xlab("") + ylab("Log2 gene count (MAD)")
print(p)

#p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=percent.mt_MAD)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_violin() + scale_fill_discrete(name='') + guides(col=FALSE)# + geom_point(position=position_jitterdodge(jitter.width=0.1), alpha=0.2)
p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=percent.mt_MAD, fill=Immune_Cell)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_boxplot() + scale_fill_discrete(name='') + guides(col=FALSE)# + geom_point(position=position_jitterdodge(jitter.width=0.1), alpha=0.2)
p = p + geom_hline(yintercept = 5, linetype=5)
p = p + geom_vline(aes(xintercept = factor(gsub("P018-PGL3", "P018-PGL1",Sample)) %>% as.integer() %>% (function(x) {x+0.5})), size=0.25)
p = p + xlab("") + ylab("Percent mitochondrial counts (MAD)")
p = p + ylim(-3,20)
print(p)

#p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=Scrublet_MAD)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_violin() + scale_fill_discrete(name='') + guides(col=FALSE)# + geom_point(position=position_jitterdodge(jitter.width=0.1), alpha=0.2)
p = ggplot(qc_metadata, aes(x=gsub("E207", "E207 (removed)", Sample), y=Scrublet_MAD, fill=Immune_Cell)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_boxplot() + scale_fill_discrete(name='') + guides(col=FALSE)# + geom_point(position=position_jitterdodge(jitter.width=0.1), alpha=0.2)
p = p + geom_hline(yintercept = 2, linetype=5)
p = p + geom_vline(aes(xintercept = factor(gsub("P018-PGL3", "P018-PGL1",Sample)) %>% as.integer() %>% (function(x) {x+0.5})), size=0.25)
p = p + xlab("") + ylab("Scrublet (MAD)")
print(p)

dev.off()


#genes = lapply(pseudobulk.normal.signatures, head, n=2) %>% unlist
#genes = pseudobulk.normal.signatures$Fibroblast3
f = order(sn.PPGLs$Cell_Type)[order(sn.PPGLs$Cell_Supertype[order(sn.PPGLs$Cell_Type)])] %>% setdiff(which(sn.PPGLs$Cell_Supertype %in% c("Tumour", "Adrenocortical cells", "Chromaffin cells", "Sustentacular cells")))
x = do.call(cbind, lapply(unique(sn.PPGLs$Cell_Type[f]), function(cell_type){
  rowSums(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Cell_Type==cell_type])
})) %>% LogNormalize
y = do.call(cbind, lapply(unique((sn.PPGLs$Cell_Type)[f]), function(cell_type){
  rowMeans(sn.PPGLs@assays$RNA@counts[,sn.PPGLs$Cell_Type==cell_type] > 0)
}))
colnames(x) = gsub("neg ", "- ", gsub("pos ", "+ ", unique((sn.PPGLs$Cell_Type)[f])))
colnames(y) = colnames(x)
z = data.frame(`Cell Type` = sn.PPGLs$Cell_Supertype[match(unique((sn.PPGLs$Cell_Type)[f]), sn.PPGLs$Cell_Type)], check.names=F)
hm.annot = HeatmapAnnotation(df = z, which="column", col = list("Cell Type" = cell.cols), show_legend = FALSE)

genes = c("FLT1",  # Endothelial cells
          "PROX1", # Lymphatic
          "DLL4", # Tip-like/Stalk-like
          "PDGFRB", # All Fibroblasts
          "ACSM3", "FMO2", "MFAP5",
          "POSTN", # Myofibroblast
          "MCAM", # Pericyte
          "ACTA2", # Smooth muscle
          "ITGAM", "FCGR3A", # Monocytes
          "HLA-DRA", "CD1C", "IDO1", # DCs
          "FCGR2A", # Macrophage
          "TPSAB1", # Mast cells
          "CD2",
          "CD4", "CD8A", "GZMB", "FCRL6", "NCAM1", "FOXP3", # T cells
          "MS4A1", "MZB1", "FCRL5" # B cells
          )

#genes = c("CD34", "PROM1")

hm = HeatmapDotPlot(colour = x[genes %>% intersect(rownames(x)),], size = y[genes %>% intersect(rownames(x)),],
                     scale = TRUE, cell.size = 0.5,
                     cluster_columns = F, cluster_rows=F, bottom_annotation=hm.annot,
                     col=colorRamp2(c(-3,0,3),c("blue","grey90","red")),
                     name="Scaled expression", row_names_gp=gpar(fontface="italic"),
                     column_split=z[,1], column_title_rot=90, norm.size = FALSE)
pdf("Raw figures/Cell_subtypes_dotplot.pdf", width=12, height=12)
do.call(draw, hm)
dev.off()

