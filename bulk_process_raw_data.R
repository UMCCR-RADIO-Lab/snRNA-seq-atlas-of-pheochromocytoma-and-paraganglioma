suppressMessages(expr = {
  library(AnnotationDbi)
  library(dplyr)
})

args = commandArgs(trailingOnly=F)

collapse_to_symbols = function(x, annotation, aggr.fun="mean"){
  
  annotation.package =  paste0(annotation, ".db")
  
  if (is.na((grep(pattern = paste0("^",as.character(annotation.package), "$"),
                  installed.packages())[1]))) {
    BiocManager::install(annotation.package)
  }
  library(annotation.package, character.only=TRUE)
  
  probe.to.symbol = toTable(toggleProbes(value = "all",
                                         x = eval(as.symbol(paste0(annotation,
                                                                   "SYMBOL")))))
  
  
  symbols = probe.to.symbol[match(rownames(x), probe.to.symbol[,1]),2]
  
  f = which(!is.na(symbols) & !(duplicated(symbols) | duplicated(symbols, fromLast=T)) & symbols != "")
  f2 = which(!is.na(symbols) & (duplicated(symbols) | duplicated(symbols, fromLast=T)) & symbols != "")
  
  # Aggregate probes to gene symbol level
  x.single = x[f,]
  rownames(x.single) = symbols[f]
  
  
  x.dup = x[f2,]
  x.dup = data.frame(Symbol = symbols[f2],
                     x.dup) %>% group_by(Symbol) %>%
    summarise_all(eval(as.symbol(aggr.fun)))
  x.dup = as.data.frame(x.dup)
  rownames(x.dup) = x.dup[,1]
  x.dup = as.matrix(x.dup[,-1])
  
  x = rbind(x.single, x.dup)
  x = data.frame(Gene=rownames(x), x, stringsAsFactors = F)
  
  return(x)
}

proc_affy = function(dir,
                     out="expression.tsv",
                     annotation="hgug133a",
                     aggr.fun = "mean"){
  library(affy)

  data.raw = ReadAffy(filenames=list.files(dir,
                                           full.names=FALSE,
                                           pattern="\\.CEL$"),
                      celfile.path=dir)
  data.rma = rma(data.raw)
  
  x = collapse_to_symbols(exprs(data.rma), annotation, aggr.fun)
  colnames(x) = gsub("\\.CEL$", "", colnames(x))
  write.table(x, file=out, quote=F, sep="\t", col.names=T, row.names=F)
  
}

proc_2c = function(dir,
                   out="expression.tsv",
                   annotation="hgug4112a",
                   aggr.fun = "mean"){
  library(limma)
  tmp = read.maimages(list.files(dir,full.names = T),source = "agilent")
  tmp = backgroundCorrect(tmp, method="normexp", offset=5)
  tmp = normalizeWithinArrays(tmp, method="loess")
  tmp = normalizeBetweenArrays(tmp)
  tmp = avereps(tmp, ID=tmp$genes$ProbeName)
  rownames(tmp$genes) = tmp$genes$ProbeName
  new_eset = new("MIAME")
  new_eset = new("ExpressionSet", exprs = as.matrix(tmp$A),
                 phenoData = new("AnnotatedDataFrame", data = tmp$targets),
                 featureData = new("AnnotatedDataFrame", data = tmp$genes))
  sampleNames(new_eset) = basename(sampleNames(new_eset))
  exprs(new_eset) = exprs(new_eset) - 2
  
  x = collapse_to_symbols(exprs(new_eset), annotation, aggr.fun)
  colnames(x) = gsub("\\.CEL$", "", colnames(x))
  write.table(x, file=out, quote=F, sep="\t", col.names=T, row.names=F)
}


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
      values = NULL
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
                 defaults = list(out="expression.tsv",
                                 type="affy",
                                 annotation="hgug4112a",
                                 aggr = "mean"))

if(!is.null(args[["dir"]])){
  dir = args[["dir"]]
  annotation = args[["annotation"]]
  type = args[["type"]]
  out = args[["out"]]
  aggr = args[["aggr"]]
  if(type == "affy"){
    proc_affy(dir, out, annotation, aggr)
  }
  if(type == "2c"){
    proc_2c(dir, out, annotation, aggr)
  }
  cat(paste0("Finished! (", out, ")\n"))
} else {
  cat("No input directory\n")
}


