library(GEOquery)




GSE2841 = getGEO(filename="data/series/GSE19987-GPL96_series_matrix.txt.gz",
                 GSEMatrix=T, AnnotGPL=F, getGPL=FALSE)
GSE19987 = getGEO(filename="data/series/GSE19987-GPL571_series_matrix.txt.gz",
                  GSEMatrix=T, AnnotGPL=F, getGPL=FALSE)
GSE67066 = getGEO(filename="data/series/GSE67066_series_matrix.txt.gz",
                  GSEMatrix=T, AnnotGPL=F, getGPL=FALSE)
GSE19422 = getGEO(filename="data/series/GSE19422_series_matrix.txt.gz",
                  GSEMatrix=T, AnnotGPL=F, getGPL=FALSE)
GSE51081 = getGEO(filename="data/series/GSE51081_series_matrix.txt.gz",
                  GSEMatrix=T, AnnotGPL=F, getGPL=FALSE)
GSE50442 = getGEO(filename="data/series/GSE50442_series_matrix.txt.gz",
                  GSEMatrix=T, AnnotGPL=F, getGPL=FALSE)

GSE39716 = getGEO(filename="data/series/GSE39716_series_matrix.txt.gz", GSEMatrix=T, AnnotGPL=F,getGPL=FALSE)
EMTAB733_annot = read.delim("data/annot/E-MTAB-733.sdrf.tsv",sep="\t",header=T,as.is=T)
EMTAB733_annot$Hybridization.Name = trimws(EMTAB733_annot$Hybridization.Name)
GSE19422_annot = read.delim("data/annot/GSE51081_readme_79re-analyzedSamples.txt",sep="\t",header=T,as.is=T)
GSE19422_annot = GSE19422_annot[!GSE19422_annot$Sample.name=="",]
flynn_annot = read.delim("data/annot/Flynn_Annotation.csv",sep=",",header=T,as.is=T)
flynn_annot = flynn_annot[!flynn_annot$TumourID == "",]
tcga_annot = read.delim("data/annot/TCGA_Annotation.csv", sep=",", header=T,as.is=T)
tcga_annot2 = read.delim("data/annot/TCGA_annotation_updated.tsv", sep="\t", header=T, as.is=T)

dahia_annot = read.delim("data/annot/Dahia-updated_annotation_with_class.tsv",sep="\t",header=T,as.is=T)
dahia_annot = dahia_annot[!duplicated(dahia_annot),]
dahia_annot2 = read.delim("data/annot/dahia_cel_file_key.txt", sep="\t", header=T, as.is=T)
dahia_annot2 = dahia_annot2[!dahia_annot2$Scan.Name=="CS2006121839",] # already in GEO as GSM499577


castro_vega_annot = read.delim("data/annot/castrovega_annot.txt", sep="\t", header=T, as.is=T, skip = 2)
castro_vega_annot = castro_vega_annot[,-grep("^X$|^X\\.",colnames(castro_vega_annot))]

annot = data.frame(ID = c(colnames(GSE2841),
                              colnames(GSE19987),
                              paste0(dahia_annot2$Scan.Name, ".CEL"),
                              colnames(GSE67066),
                              colnames(GSE39716),
                              EMTAB733_annot$Hybridization.Name,
                              colnames(GSE19422),
                              colnames(GSE51081),
                              flynn_annot$TumourID,
                              tcga_annot2$CASE_ID,
                              colnames(GSE50442)
                              ),
                   Dataset = c(
                     rep("GSE2841", ncol(GSE2841)),
                     rep("GSE19987", ncol(GSE19987)),
                     rep("GSE19987 (missing)", nrow(dahia_annot2)),
                     rep("GSE67066", ncol(GSE67066)),
                     rep("GSE39716", ncol(GSE39716)),
                     rep("E-MTAB-733", nrow(EMTAB733_annot)),
                     rep("GSE19422", ncol(GSE19422)),
                     rep("GSE51081", ncol(GSE51081)),
                     rep("Flynn", nrow(flynn_annot)),
                     rep("TCGA", nrow(tcga_annot2)),
                     rep("GSE50442", ncol(GSE50442))
                               ), stringsAsFactors = F)
#annot$ID = trimws(annot$ID)
annot[,c("Sample", "Accession_ID", "PCPG", "Genotype_RAW", "Location_RAW", "Malignancy_RAW", "Sex", "Age", "Adrenal_Activity", "Purity")] = NA
rownames(annot) = annot$ID

# GSE2841
f = match(colnames(GSE2841), annot$ID)
annot$Sample[f] = gsub("\ .*$", "", as.character(pData(GSE2841)$title))
annot$Accession_ID[f] = as.character(pData(GSE2841)$geo_accession)
annot$PCPG[f] = "PCC"
annot$Genotype_RAW[f] = as.character(pData(GSE2841)$characteristics_ch1.1)
annot$Location_RAW[f] = as.character(pData(GSE2841)$characteristics_ch1.2)

# GSE19987
f = match(colnames(GSE19987), annot$ID)
annot$Sample[f] = gsub("\ .*$", "", as.character(pData(GSE19987)$title))
annot$Accession_ID[f] = as.character(pData(GSE19987)$geo_accession)
annot$PCPG[f] = "PCC"
annot$Genotype_RAW[f] = as.character(pData(GSE19987)$characteristics_ch1.1)
annot$Location_RAW[f] = as.character(pData(GSE19987)$characteristics_ch1.2)



# Dahia missing
annot$Sample[match(paste0(dahia_annot2$Scan.Name, ".CEL"), annot$ID)] = paste0("P", gsub("^.* ", "", dahia_annot2$Sample))
dahia_annot = dahia_annot[dahia_annot$sample.ID %in% annot$Sample,]
f = match(dahia_annot$sample.ID, annot$Sample)
annot$Genotype_RAW[f] = dahia_annot$Gene




# GSE19422
f = match(colnames(GSE19422), annot$ID)
annot$Sample[f] = as.character(pData(GSE19422)$description)
annot$Accession_ID[f] = colnames(GSE19422)
annot$PCPG[f] = gsub("\ .*$", "", as.character(pData(GSE19422)$title))
annot$Genotype_RAW[f] = gsub("_.*$", "", GSE19422_annot$Sample.name[match(annot$Sample[f], gsub("-.*$", "",GSE19422_annot$description.1))])

# GSE51081
f = match(colnames(GSE51081), annot$ID)
annot$Sample[f] = gsub("^.*(|).*$", "", as.character(GSE51081$title))
annot$Accession_ID[f] =  colnames(GSE51081)
annot$PCPG[f] = "PCC"
annot$Genotype_RAW[f] = as.character(GSE51081$`mutation:ch1`)

# E-MTAB-733

f = match(EMTAB733_annot$Hybridization.Name, annot$ID)
annot$Sample[f] = EMTAB733_annot$Source.Name
annot$Accession_ID[f] = EMTAB733_annot$Protocol.REF.3
#annot$PCPG[f] = NA
annot$Sex[f] = EMTAB733_annot$Characteristics..Sex.
annot$Genotype_RAW[f] = EMTAB733_annot$Characteristics..Genotype.

sample.numbers = as.numeric(gsub("^.*_", "", EMTAB733_annot$Source.Name))
castro_vega_annot = castro_vega_annot[castro_vega_annot$Sample.. %in% sample.numbers,]
f = match(castro_vega_annot$Sample.., ifelse(annot$Dataset == "E-MTAB-733",
                                             as.numeric(gsub("^.*_", "", annot$Sample)), NA))
annot$Malignancy_RAW[f] = castro_vega_annot$Benign..Malignant
annot$Location_RAW[f] = castro_vega_annot$Tumor.Location
annot$Sex[f] = castro_vega_annot$Sex
annot$Age[f] = castro_vega_annot$Age.at.diagnosis..years.
annot$Genotype_RAW[f] = paste(castro_vega_annot$Germline.Mutation, castro_vega_annot$Somatic.Mutation, sep=";")


# GSE39716
f = match(colnames(GSE39716), annot$ID)
annot$Sample[f] = as.character(GSE39716$title)
annot$Accession_ID[f] = colnames(GSE39716)
annot$PCPG[f] = as.character(GSE39716$source_name_ch1)
annot$Age[f] = as.character(GSE39716$`age(yrs):ch1`)
annot$Sex[f] = as.character(GSE39716$`gender:ch1`)
annot$Location_RAW[f] = as.character(GSE39716$characteristics_ch1.4)
annot$Genotype_RAW[f] = as.character(GSE39716$`mutation:ch1`)
annot$Malignancy_RAW[f] = as.character(GSE39716$characteristics_ch1.4)



# GSE67066
f = match(colnames(GSE67066), annot$ID)
annot$Sample[f] = as.character(GSE67066$title)
annot$Accession_ID[f] = colnames(GSE67066)
annot$PCPG[f] = as.character(GSE67066$source_name_ch1)
annot$Malignancy_RAW[f] = as.character(GSE67066$characteristics_ch1.1)
annot$Genotype_RAW[f] = as.character(GSE67066$`mutation type:ch1`)

# Flynn
f = match(flynn_annot$TumourID, annot$ID)
annot$Sample[f] = flynn_annot$TumourID
annot$Purity[f] = flynn_annot$Tumour.purity.ASCAT.
annot$Genotype_RAW[f] = flynn_annot$Driver.Gene
annot$Location_RAW[f] = flynn_annot$Diagnosis
annot$Malignancy_RAW[f] = flynn_annot$MalignancyClass
annot$PCPG[f] = flynn_annot$Type
annot$Age[f] = flynn_annot$Age.at.Surgery
annot$Sex[f] = flynn_annot$Sex
annot$Adrenal_Activity[f] = flynn_annot$Presence.of.elevated.catecholamines



# TCGA
f = match(tcga_annot2$CASE_ID, annot$ID)
annot$Sample[f] = tcga_annot2$CASE_ID
annot$Age[f] = tcga_annot2$AGE
annot$PCPG[f] = tcga_annot2$CANCER_TYPE_DETAILED

annot$Sex[f] = tcga_annot2$SEX
f2 = match(tcga_annot2$CASE_ID, substr(tcga_annot$Sample.ID, 1, 15))
annot$Genotype_RAW[f] = tcga_annot$Gene.Driver[f2]
annot$Purity[f] = tcga_annot$Purity[f2]
cols = grep("\\.Secreting", colnames(tcga_annot), value=T)
annot$Adrenal_Activity[f] = apply(tcga_annot[f2,cols], 1, function(x){
  paste(cols[which(x == "Yes")], collapse=",")
})
annot$Malignancy_RAW[f] = tcga_annot$Metastatic.Disease[f2]
annot$Location_RAW[f] = tcga_annot$Tumor.Location[f2]
annot$Location_RAW[f] = ifelse(is.na(annot$Location_RAW[f]), tcga_annot2$TUMOR_SITE, annot$Location_RAW[f])

# GSE50442
f = match(colnames(GSE50442), annot$ID)
annot$Sample[f] = as.character(GSE50442$title)
annot$Accession_ID[f] = colnames(GSE50442)
annot$Location_RAW[f] = as.character(GSE50442$`tissue:ch1`)
annot$PCPG[f] = as.character(GSE50442$`sample type:ch1`)
annot$Malignancy_RAW[f] = as.character(GSE50442$characteristics_ch1.2)
annot$Age[f] = as.character(GSE50442$`age at surgery(yrs):ch1`)
annot$Sex[f] = as.character(GSE50442$`gender:ch1`)
annot$Genotype_RAW[f] = as.character(GSE50442$description)





write.table(annot, file="data/annot/merged_pdata.tsv", sep="\t", quote=F, col.names=T, row.names=F)


