library(preprocessCore)
library(parallel)
library(sva)
library(limma)
## Collect arguments
args <- commandArgs(trailingOnly = TRUE)
#args = unlist(strsplit(c("--microarray", list.files("data/microarrays/", full.names=T, pattern=".tsv"),
#             "--rnaseq", list.files("data/RNASeq/", full.names=T, pattern=".tsv"),
#             "--out", "data/merged_unfiltered.tsv"), " "))

## Help section
help = '
      Merge datasets and remove batch effects

Arguments:
--rnaseq rnaseq.tsv - tab separated rnaseq files (gene name in first column)
                      can be multiple files
--microarray array.tsv - tab separated microarray files (gene name in first column)
                              can be multiple files
--out merged.tsv    - output file
--filter badsamples.txt - file with samples to be filtered out, no column names
--factors genotypes.txt - file with factors to be preserved
                          first column should be sample names (as read by R), need not be sorted
                          other columns can be factors used to create design matrix
                          can be multiple files with same column names
                          missing values can work, they will be given 0 weight in lmFit
--quantilenormfirst - quantile normalise all batches only once, regardless of missing genes
--help              - print this text

Example:
./Rscript merge_datasets.R -microarray data/E-MTAB-733.tsv -rnaseq data/TCGA.tsv -out data/merged.tsv \n\n'


if("--help" %in% args | "-h" %in% args | length(args) == 0) {
  cat(help)
  
  q(save="no")
}

## Parse arguments
parseArgs <- function(x, defaults = list()) {
  args = defaults
  x = unlist(strsplit(x, "="))
  is.key = grepl("^-", x)
  key = ""
  for(i in 1:length(x)){
    arg = x[i]
    if(is.key[i]){
      key = gsub("^-*", "", arg)
      if(i == length(x) | is.key[min(i+1, length(x))]){
        args[[key]] = TRUE
      }
      values = args[[key]]
    } else {
      values = c(values, arg)
      if(i == length(x) | is.key[min(i+1, length(x))]){
        args[[key]] = values
      }
    }
  }
  
  return(args)
}

args = parseArgs(args,
                 defaults = list("quantilenormfirst" = FALSE,
                                 "saveintermediate" = FALSE))

# Samples to remove before within-batch processing
filtered.samples = args[["filter"]]
if(!is.null(filtered.samples[1])){
  #filtered.samples = read.delim(filtered.samples, sep="\t", header=F, as.is=T)[,1]
  filtered.samples =  do.call(c, lapply(filtered.samples, function(x){
    read.delim(x, sep="\t", header=F, as.is=T)[,1]
  }))
}
cat("Filtering:\n")
print(filtered.samples)

output.file = c(args[["o"]], args[["out"]], args[["output"]])[1]
if(is.null(output.file)){
  stop("")
}

# Process microarray datasets
microarrays = mclapply(args[["microarray"]], function(x){
  cat(paste0("Loading ", x, "\n"))
  data = read.delim(x, sep="\t", header=T, as.is=T, row.names=1)
  return(data[,setdiff(colnames(data), filtered.samples)])
}, mc.cores = 4)
names(microarrays) = gsub("\\.txt|\\.tsv", "", basename(args[["microarray"]]))

# Process RNASeq datasets
rnaseq = mclapply(args[["rnaseq"]], function(x){
  cat(paste0("Loading ", x, "\n"))
  data = read.delim(x, sep="\t", header=T, as.is=T, row.names=1)
  return(data[,setdiff(colnames(data), filtered.samples)])
}, mc.cores=4)
names(rnaseq) = gsub("\\.txt|\\.tsv", "", basename(args[["rnaseq"]]))

# Read factors
# 1st column of file = Sample IDs
factors = args[["factors"]]
if(!is.null(factors[1])){
  factors =  do.call(rbind, lapply(factors, function(x){
    read.delim(x, sep="\t", header=T, as.is=T, row.names=1)
  }))
}

normalize.quantiles.keepdimnames = function(x, func = normalize.quantiles, ...){
  x2 = func(x, ...)
  dimnames(x2) = dimnames(x)
  return(x2)
}

removeBatchEffectModel = function (x, batch = NULL, batch2 = NULL, covariates = NULL, 
          design = matrix(1, ncol(x), 1), ...) 
{
  if (is.null(batch) && is.null(batch2) && is.null(covariates)) 
    return(as.matrix(x))
  if (!is.null(batch)) {
    batch <- as.factor(batch)
    contrasts(batch) <- contr.sum(levels(batch))
    batch <- model.matrix(~batch)[, -1, drop = FALSE]
  }
  if (!is.null(batch2)) {
    batch2 <- as.factor(batch2)
    contrasts(batch2) <- contr.sum(levels(batch2))
    batch2 <- model.matrix(~batch2)[, -1, drop = FALSE]
  }
  if (!is.null(covariates)) 
    covariates <- as.matrix(covariates)
  X.batch <- cbind(batch, batch2, covariates)
  fit <- lmFit(x, cbind(design, X.batch), ...)
  fit = (as.data.frame(fit$coefficients))
  rownames(fit) = rownames(x)
  return(fit)
}

mergeExpression = function(x = list()){
  cols = sort(unique(unlist(lapply(x, colnames))))
  rows = sort(unique(unlist(lapply(x, rownames))))
  #m = data.frame(NA, row.names = rows, col.names=cols)
  m = matrix(NA, nrow = length(rows), ncol=length(cols), dimnames = list(rows, cols))
  for(y in x){
    m[rownames(y),colnames(y)] = as.matrix(y)
  }
  return(m)
}

# Quantile-normalise, then remove batch effect using ComBat
renormBatch = function(x, batches, genes.keep=NULL,
                        x2 = NULL, batches2=NULL,
                        factors=list(),
                       quantile.norm=TRUE, return.coefficients=FALSE){
  
  if(quantile.norm){
    #x.n = normalize.quantiles(x, copy=T)
    #dimnames(x.n) = dimnames(x) # annoyingly normalize.quantiles wipes dimnames
    #x = x.n; rm(x.n)
    x = normalize.quantiles.keepdimnames(x, copy=TRUE)
    
    if(!is.null(dim(x2))){
      if(ncol(x2) > 0){
        #x2.n = normalize.quantiles.use.target(x2, target=x[,1], copy = TRUE)
        #dimnames(x2.n) = dimnames(x2)
        #x2 = x2.n; rm(x2.n)
        x2 = normalize.quantiles.keepdimnames(x, func=normalize.quantiles.use.target,
                                              target=x[,1], copy=TRUE)
        
        x = cbind(x, x2)
        batches = c(batches, batches2)
      }
    }
  } else {
    x = cbind(x, x2)
    batches = c(batches, batches2)
  }
  
  
  if(!is.null(genes.keep[1])){
    x = x[genes.keep,,drop=F]
  }
  
  if(length(factors)>0){
    factors = factors[match(colnames(x), rownames(factors)),,drop=F]
    if(sum(is.na(factors)) > 0){
      weights = ifelse(apply(is.na(factors), 1, sum) > 0, 0, 1)
      factors[is.na(factors)] = 0 # dummy values
      mod = model.matrix(as.formula(paste("~0+",paste(colnames(factors), collapse = "+"))), data = factors)
      
      if(return.coefficients){
        x.norm = removeBatchEffectModel(x, batch=batches, design=mod, weights=weights)
      } else {
        x.norm = removeBatchEffect(x, batch=batches, design=mod, weights=weights)
      }
      
      
    } else {
      mod = model.matrix(as.formula(paste("~0+",paste(colnames(factors), collapse = " + "))), data = factors)
      if(return.coefficients){
        x.norm = removeBatchEffectModel(x, batch=batches, design=mod)
      } else {
        x.norm = removeBatchEffect(x, batch=batches, design=mod)
      }
    }
    
  } else {
    mod = NULL
    if(return.coefficients){
      x.norm = removeBatchEffectModel(x, batch=batches)
    } else {
      x.norm = removeBatchEffect(x, batch=batches)
    }
  }
  #x.norm = ComBat(dat=x, batch=batches, mod=mod, par.prior=TRUE, prior.plots=FALSE)
  
  return(x.norm)
}

mergeAll = function(microarrays, rnaseq, factors = NULL, quantilenormfirst=FALSE,
                    saveintermediate=FALSE, savemodel=FALSE){
  
  microarray.genes = unique(unlist(lapply(microarrays, rownames)))
  # Ignore genes only found in RNA-Seq
  # to avoid any quantile normalisation issues
  rnaseq = lapply(rnaseq, function(x){
    x[intersect(rownames(x), microarray.genes),]
  })
  
  microarrays.batches = unlist(lapply(names(microarrays), function(batch){
    setNames(rep(batch, ncol(microarrays[[batch]])),
             colnames(microarrays[[batch]]))
  }))
  rnaseq.batches = unlist(lapply(names(rnaseq), function(batch){
    setNames(rep(batch, ncol(rnaseq[[batch]])),
             colnames(rnaseq[[batch]]))
  }))
  
  all.merged = mergeExpression(c(microarrays, rnaseq))
  all.batches = c(rnaseq.batches, microarrays.batches)[colnames(all.merged)]
  is.rnaseq = colnames(all.merged) %in% unlist(lapply(rnaseq, colnames))
  
  if(quantilenormfirst){
    #microarrays.quantilenormed = normalize.quantiles(all.merged[,!is.rnaseq], copy=TRUE)
    microarrays.quantilenormed = normalize.quantiles.keepdimnames(all.merged[,!is.rnaseq], copy=TRUE)
    rnaseq.quantilenormed = normalize.quantiles.keepdimnames(all.merged[,is.rnaseq],
                                                             func = normalize.quantiles.use.target,
                                                             target=normalize.quantiles.determine.target(all.merged[,!is.rnaseq]),
                                                             copy=TRUE)
    all.merged = cbind(microarrays.quantilenormed, rnaseq.quantilenormed)[,colnames(all.merged)]
    rm(microarrays.quantilenormed,
       rnaseq.quantilenormed)
    gc()
    if(!is.null(qn.output.file)){
      write.table(all.merged, file=qn.output.file, sep="\t", quote=F,
                  row.names=TRUE, col.names=NA)
      write.table(data.frame(Sample = colnames(all.merged), Batch = all.batches),
                  file = gsub("\\.tsv", "_batches.tsv", qn.output.file), sep="\t",
                  quote=F, row.names=F, col.names=T)
    }
    
  }
  
  #snowparam <- SnowParam(workers = 4, type = "SOCK")
  #register(snowparam, default = TRUE)
  
  lm.coefficients = NULL
  
  cat("Normalising...\n")
  combinations = unique(!is.na(all.merged)) # all gene/batch combinations
  for(i in 1:nrow(combinations)){
    cat(paste0("\r", round(100*((i-1)/nrow(combinations)), digits=1), "%    "))
    sink("/dev/null")
    suppressWarnings({
      suppressMessages({
        
        sample.subset = combinations[i,]
        
        subset.merged = all.merged[,sample.subset]
        subset.merged = subset.merged[apply(is.na(subset.merged), 1, sum) == 0,]
        subset.batches = all.batches[sample.subset]
        
        # only save values that are *only* not missing within this subset
        # But normalise against all genes not missing within those samples
        genes.keep = setdiff(rownames(subset.merged),
                             rownames(all.merged)[apply(!is.na(all.merged[,!sample.subset]), 1, sum) > 0])
        if(length(unique(subset.batches)) > 1 & length(genes.keep) > 0){
          subset.is.rnaseq = is.rnaseq[sample.subset]
          subset.normed = renormBatch(x = subset.merged[,!subset.is.rnaseq],
                                       batches = subset.batches[!subset.is.rnaseq],
                                       x2 = subset.merged[,subset.is.rnaseq],
                                       batches2 = subset.batches[subset.is.rnaseq],
                                       factors = factors,
                                       genes.keep = genes.keep,
                                      quantile.norm = !quantilenormfirst)
          
          all.merged[genes.keep, colnames(subset.normed)] = subset.normed[genes.keep,]
          if(savemodel){
             a = rownames(lm.coefficients)
             lm.coefficients = plyr::rbind.fill(lm.coefficients,
                                    renormBatch(x = subset.merged[,!subset.is.rnaseq],
                                                batches = subset.batches[!subset.is.rnaseq],
                                                x2 = subset.merged[,subset.is.rnaseq],
                                                batches2 = subset.batches[subset.is.rnaseq],
                                                factors = factors,
                                                genes.keep = genes.keep,
                                                quantile.norm = !quantilenormfirst, return.coefficients =TRUE)[genes.keep,])
            rownames(lm.coefficients)=c(a, genes.keep)
          }
        }
      })
    })
    sink()
    
    
    
    
  }
  cat("\r100.0%   \n")
  cat("Finished normalising\n")
  if(savemodel){write.table(lm.coefficients, file=model.output.file, sep='\t', quote=F, col.names=NA, row.names=T)}
  return(all.merged)
}

quantilenormfirst = args[["quantilenormfirst"]]
saveintermediate = args[["saveintermediate"]]
if(saveintermediate){
  output.file = gsub("\\....$|\\...$", "", output.file)
  qn.output.file = paste0(output.file, "_quantilenormed.tsv")
  output.file = paste0(output.file, "_batchnormed.tsv")
} else {
  qn.output.file = NULL
}
savemodel = args[["savemodel"]]
if(savemodel){model.output.file = paste0(gsub("\\.tsv", "", output.file), ".coef.tsv"); print(model.output.file)}

merged = mergeAll(microarrays, rnaseq, factors, quantilenormfirst, qn.output.file, savemodel = savemodel)

cat(paste0("Writing ", output.file, "\n"))
write.table(data.frame(Gene = rownames(merged), merged),
            file=output.file, quote=F, sep="\t",
            row.names=F, col.names=T)
cat("Finished!\n")
