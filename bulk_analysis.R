library(tidyverse)
library(limma)
library(edgeR)

# Change to the google drive dir
setwd("~/Google Drive/PPGL single cell paper/")

bulk.mat = read.delim("Raw tables/Bulk/bulk_expression_quantile_normalised.tsv", sep='\t', header=T, as.is=T, check.names=F, row.names=1)
bulk.mat2 = read.delim("Raw tables/Bulk/bulk_expression_batch_normalised.tsv", sep='\t', header=T, as.is=T, check.names=F, row.names=1)

bulk.meta = read.delim("Raw tables/Bulk/bulk_metadata.tsv", sep='\t', header=T, as.is=T)
bulk.meta$Genotype = factor(bulk.meta$Genotype,
                            levels=c("Normal",
                                     "RET",
                                     "HRAS",
                                     "NF1",
                                     "BRAF",
                                     "TMEM127",
                                     "NGFR",
                                     "MAX",
                                     "CSDE1",
                                     "H3F3A",
                                     "KIF1B",
                                     "IDH1",
                                     "MAML3",
                                     "FH",
                                     "EPAS1",
                                     "EGLN1",
                                     "VHL",
                                     "SDHA",
                                     "SDHB",
                                     "SDHC",
                                     "SDHD",
                                     "SETD2",
                                     "Unknown"))

# UMAP/Consensus Clustering

bulk.hvg = bulk.mat2 %>% (function(x, nmad=3){
  v = apply(x, 1, function(x){
    sd(x, na.rm=T)
  })
  f = apply(is.na(x), 1, sum) == 0
  v.t = median(v[f]) + nmad*mad(v[f])
  f = f & v >= v.t
  x=t(x[f,])
  return(x)
})

run_umap = function(x, do.scale=FALSE, ...){
  if(do.scale){
    centers = apply(x,2,mean, na.rm=T)
    x=scale(x, scale=FALSE, center=TRUE)
  }
  u = uwot::umap(x, ...)
  #u = uwot::umap(as.dist(1-cor(x)),  metric = "cosine", nn_method="annoy", n_epochs = n_epochs, ...)
  u$features=colnames(x)
  if(do.scale){
    u$centers = centers
  }
  return(u)
}


cosine.dist = function(x,y){
  f = !is.na(x) & !is.na(y)
  x = x[f]
  y=y[f]
  a=sum(x*y)
  b=sum(x**2)
  c=sum(y**2)
  if(b==0 & c ==0){return(0)}
  if(b == 0 | c == 0){return(1)}
  return(1 - a/sqrt(b*c))
}
cosine.knn.impute = function(x, y=NULL, k=15, aggr.fun = mean){
  xo = x
  if(is.null(dim(y))){
    y = x
  } else {
    y = y[match(rownames(x), rownames(y)),]
  }
  
  f = which(colSums(is.na(x)) > 0)
  for(j in f){
    cat(paste0(as.integer(100*match(j,f)/length(f)), "%\r"))
    
    f2 = which(is.na(x[,j]))
    if(length(f2) > 0){
      nn = order(apply(y,2,cosine.dist, y=x[,j]), decreasing = T)
    }
    xo[f2,j] = sapply(f2, function(i){
      aggr.fun(head(na.omit(as.numeric(y[i,nn])), k), na.rm=T)
    })
  }
  return(xo)
}

set.seed(1);bulk.umap.model = run_umap(bulk.hvg, do.scale=TRUE,
                                       ret_model=TRUE, min_dist=0.1,
                                       n_epochs = 1000, metric='cosine',
                                       nn_method='annoy', n_neighbors=15)

bulk.dist = bulk.hvg %>% (function(x) {
  x=x %>% scale(scale=F,center=T) %>% t
  as.dist(do.call(cbind, lapply(1:ncol(x), function(i){
    sapply(1:ncol(x), function(j){cosine.dist(x[,i], x[,j])}) %>% setNames(colnames(x))
  }) %>% setNames(colnames(x))))
})
bulk.ccp = ConsensusClusterPlus(d = as.dist(bulk.dist), distance='pearson', reps = 1000,
                                maxK = 12,pItem = 0.7, clusterAlg = 'hc',
                                title = 'Raw tables/Bulk/Bulk_CCP',
                                innerLinkage = 'ward.D2', finalLinkage="ward.D2",
                                seed = 1, plot = 'pdf',
                                corUse = 'pairwise.complete.obs')

bulk.consensus = calcICL(bulk.ccp)
bulk.meta$Consensus = bulk.consensus$itemConsensus$itemConsensus[match(paste(colnames(bulk.mat2), bulk.meta$ConsensusCluster, "9"), paste(bulk.consensus$itemConsensus$item, bulk.consensus$itemConsensus$cluster, bulk.consensus$itemConsensus$k))]
bulk.meta$ConsensusCluster = bulk.ccp[[9]]$consensusClass
bulk.meta$Cluster = c("Kinase", "Kinase", "VHL", "IDC", "PH-NOS", "Cortical admixture", "SDHx (H&N)", "SDHx", "MAML3")[bulk.meta$ConsensusCluster]

bulk.meta$Cluster[which(bulk.meta$Genotype=="Normal")] = "Normal"
bulk.meta$Cluster = factor(bulk.meta$Cluster, levels = c("Normal", "Cortical admixture", "Kinase", "IDC", "MAML3", "PH-NOS", "VHL", "SDHx", "SDHx (H&N)"))
colnames(bulk.mat) = bulk.meta$Sample
colnames(bulk.mat2) = bulk.meta$Sample

bulk.meta[,c("UMAP_1", "UMAP_2")] = bulk.umap.model$embedding

write.table(bulk.meta, file="Raw tables/Bulk/bulk_metadata.tsv", sep='\t', quote=F, col.names=T, row.names=F)

# Differential expression

# Performs differential expression between two groups
# Accounts for NA values in the coefficients
# (because bulk data has lots of missing values)
compare.groups = function(fit, groups1, groups2, trend=FALSE){
  groups1 = intersect(groups1, colnames(fit$coefficients))
  groups2 = intersect(groups2, colnames(fit$coefficients))
  fit = fit[,c(groups1, groups2)]
  all.unique.groups = !is.na(fit$coefficients) %>% as.data.frame
  unique.groups = all.unique.groups[!duplicated(all.unique.groups),,drop=F]
  
  y = do.call(rbind, lapply(1:nrow(unique.groups), function(x){
    x = unique.groups[x,]
    fit = fit[,colnames(fit$coefficients)[x]]
    groups1 = intersect(groups1, colnames(fit$coefficients))
    groups2 = intersect(groups2, colnames(fit$coefficients))
    if(length(groups1) == 0 | length(groups2) == 0){return(NULL)}
    coef = paste0("(", paste(groups1, collapse = " + "),")/", length(groups1)," - (", paste(groups2, collapse = " + "),")/", length(groups2))
    y = contrasts.fit(fit, makeContrasts(contrasts = coef, levels=colnames(fit$coefficients))) %>% eBayes(trend=trend) %>% topTable(coef=coef, sort.by = 'p', number = Inf, p.value = 1)
    y = y[rownames(y) %in% rownames(fit$coefficients)[apply(all.unique.groups, 1, function(i){sum(i!=x) == 0})],]
    return(y)
  }))
  y$gene = rownames(y)
  y = y[order(y$P.Value),]
  y$adj.P.Val = p.adjust(y$P.Value, method="BH")
  return(y)
}

bulk.design = model.matrix(~0 + Genotype + Batch, data=bulk.meta %>% filter(!Cluster == "Cortical admixture"))
colnames(bulk.design)=gsub("^Batch|^Genotype", "", colnames(bulk.design))
rownames(bulk.design)=bulk.meta %>% filter(!Cluster=="Cortical admixture") %>% pull(Sample)
bulk.design = bulk.design[,!grepl("Unknown", colnames(bulk.design))]

bulk.fit.genotype = lmFit(bulk.mat[,rownames(bulk.design)], design = bulk.design)

bulk.design.subtype = model.matrix(~0 + Cluster + Batch , data=bulk.meta)
colnames(bulk.design.subtype)=gsub("^Batch|^Cluster", "", colnames(bulk.design.subtype))
rownames(bulk.design.subtype)=bulk.meta %>% pull(Sample)
bulk.design.subtype = bulk.design.subtype[,!grepl("Unknown", colnames(bulk.design.subtype))]
colnames(bulk.design.subtype)[colnames(bulk.design.subtype) == "SDHx (H&N)"] = "SDHx_HN"
colnames(bulk.design.subtype)[colnames(bulk.design.subtype) == "Cortical admixture"] = "Cortical.admixture"
colnames(bulk.design.subtype)[colnames(bulk.design.subtype) == "PH-NOS"] = "PH_NOS"
bulk.fit.subtype = lmFit(bulk.mat, design = bulk.design.subtype)

groups = setdiff(colnames(bulk.fit.genotype$coefficients), unique(bulk.meta$Batch))
bulk.de.by.genotype = do.call(rbind, lapply(groups, function(g1){
  g2 = setdiff(groups, g1)
  y = compare.groups(bulk.fit.genotype, g1, g2, trend=TRUE)
  y$Cluster = g1
  return(y)
}))

groups = setdiff(colnames(bulk.fit.subtype$coefficients), unique(bulk.meta$Batch))
# Remove normals from the 1 vs rest
groups <- groups[-1]
bulk.de.by.subtype = do.call(rbind, lapply(groups, function(g1){
  g2 = setdiff(groups, g1)
  y = compare.groups(bulk.fit.subtype, g1, g2, trend=TRUE)
  y$Cluster = g1
  return(y)
}))

# Save
write.table(bulk.de.by.genotype, file="Raw tables/Bulk/bulk_de_by_genotype_reversed.tsv", sep='\t', quote=F, col.names=T, row.names=F)
write.table(bulk.de.by.subtype, file="Raw tables/Bulk/bulk_de_by_subtype_reversed.tsv", sep='\t', quote=F, col.names=T, row.names=F)

# Subtype-cluster association

x = table(bulk.meta$Cluster, bulk.meta$Genotype)[-1,-1]
y = do.call(cbind, lapply(1:ncol(x), function(j){
  sapply(1:nrow(x), function(i){
    m = matrix(c(sum(x[-i,]) - sum(x[-i,j]), sum(x[-i,j]),sum(x[i,]) - x[i,j], x[i,j]), nrow=2)
    t = fisher.test(m, alternative = "greater")
    return(t$p.value)
  })
}))
dimnames(y) = dimnames(x)
z = reshape::melt(y)
colnames(z) = c("Cluster", "Genotype", "p.value")
z$Cluster=as.character(z$Cluster)
z$Genotype=as.character(z$Genotype)
z$OldCluster = all.clusters[match(z$Cluster, all.clusters2)]
z$adj.p.value = p.adjust(z$p.value, method="BH")
z = cbind(z, do.call(rbind, lapply(1:nrow(z), function(i){
  data.frame(Genotype.in.cluster = x[z[i,1],z[i,2]],
             Genotype.not.in.cluster = sum(x[-match(z[i,1],rownames(x)),z[i,2]]),
             Other.in.cluster = sum(x[z[i,1],-match(z[i,2], colnames(x))]),
             Other.not.in.cluster = sum(x[-match(z[i,1],rownames(x)),-match(z[i,2], colnames(x))]))
})))
z = cbind(z, do.call(rbind, lapply(1:nrow(z), function(i){
  m = matrix(as.numeric(z[i,c(9,8,7,6)]), nrow=2)
  t = fisher.test(m, alternative="greater")
  return(data.frame(OR=t$estimate, conf.int.lower = t$conf.int[1], conf.int.upper=t$conf.int[2], stringsAsFactors = F))
})))
write.table(z, file="Raw tables/Bulk/cluster_genotype_association_p-values.tsv", sep='\t', col.names=T, row.names=F, quote=F)

# Run bulk DE by subtype comparing T vs N
bulk.design.subtype = model.matrix(~0 + Cluster + Batch , data=bulk.meta)
colnames(bulk.design.subtype)=gsub("^Batch|^Cluster", "", colnames(bulk.design.subtype))
rownames(bulk.design.subtype)=bulk.meta %>% pull(Sample)
bulk.design.subtype = bulk.design.subtype[,!grepl("Unknown", colnames(bulk.design.subtype))]
colnames(bulk.design.subtype)[colnames(bulk.design.subtype) == "SDHx (H&N)"] = "SDHx_HN"
colnames(bulk.design.subtype)[colnames(bulk.design.subtype) == "PH-NOS"] = "PH_NOS"
colnames(bulk.design.subtype)[colnames(bulk.design.subtype) == "Cortical admixture"] = "Cortical.admixture"
bulk.fit.subtype = lmFit(bulk.mat, design = bulk.design.subtype)

# Run Magnus' function comapring T vs N
groups = setdiff(colnames(bulk.fit.subtype$coefficients), unique(bulk.meta$Batch))
groups <- groups[-1]
bulk.de.by.subtype.vs.normal = do.call(rbind, lapply(groups, function(g1){
  g2 = "Normal"
  y = compare.groups(bulk.fit.subtype, g1, g2, trend=TRUE)
  y$Cluster = g1
  return(y)
}))

# Save
write.table(bulk.de.by.subtype.vs.normal, file="Raw tables/Bulk/bulk_de_by_subtpye_vs_normal.tsv", sep='\t', quote=F, col.names=T, row.names=F)

