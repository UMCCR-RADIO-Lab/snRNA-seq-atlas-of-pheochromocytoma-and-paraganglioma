library(RColorBrewer)
library(ggplot2)

# This script sets up consistent colour scheme for the pheo groups/cell types

all.genotypes = c("HRAS", "NF1","EPAS1","SDHD","MAX","FH","RET","Normal","SDHB","MAML3","TMEM127","H3F3A","VHL","Unknown","SDHA","SDHC","EGLN1","SETD2","KIF1B","MAML","NGFR","BRAF","CSDE1","IDH1")
all.cell.types = c("Tumour","Chromaffin cells",  "Adrenocortical cells", "Endothelial cells", "Fibroblasts", "Sustentacular cells", "Myeloid cells", "T/NK cells", "B cells", "Mast cells")
all.clusters = c("Kinase", "IDC", "MAML3", "PH-NOS", "VHL", "SDHx", "SDHx (H&N)", "Cortical admixture", "Normal")
all.clusters2 = c("C2A", "C2Bi", "C2Bii", "C1Bii", "C1Bi", "C1Ai", "C1Aii", "C2C", "Normal")


# Stolen from https://stackoverflow.com/questions/9563711/r-color-palettes-for-many-data-classes
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


genotype.cols = setNames(c25[-25], all.genotypes)
genotype.cols["Unknown"] = alpha('grey50', 0.25) # fade out missing genotypes
genotype.cols = genotype.cols[c(setdiff(names(genotype.cols), "Unknown"), "Unknown")]
genotype.cols["VHL"] = "darkorchid1" # remove bright yellow for white backgrounds
genotype.cols["H3F3A"] = "grey10"


cell.cols = setNames(brewer.pal(n = length(all.cell.types), name = "Set3")[c(9,10, 1:8)], all.cell.types) # reorder to make Tumour grey
cell.cols[c("Endothelial cells")] = "yellow3" # replace bright yellow for white backgrounds
cell.cols[make.names(names(cell.cols))] = unname(cell.cols) # include fixed R column names

cluster.cols = setNames(c(brewer.pal(n=8, "Set1"), "lightblue"), all.clusters)
cluster.cols["Cortical Admixture"] = cluster.cols["Cortical admixture"]
cluster.cols["SDHx_HN"] = cluster.cols["SDHx (H&N)"]
cluster.cols["SDHx"] = "yellow3" # replace bright yellow for white backgrounds
cluster.cols[make.names(names(cluster.cols))] = cluster.cols # include fixed R column names
#cluster.cols["MAML3"] = cluster.cols["MAML"]

cluster.cols[all.clusters2] = cluster.cols[all.clusters]

# ggplot scales

genotype.scale = scale_color_manual(values=genotype.cols, name="Genotype")
genotype.scale2 = scale_fill_manual(values=genotype.cols, name="Genotype")

cell.scale = scale_color_manual(values=cell.cols, name="Cell Type")
cell.scale2 = scale_fill_manual(values=cell.cols, name="Cell Type")

cluster.scale = scale_color_manual(values=cluster.cols, name="")
cluster.scale2 = scale_fill_manual(values=cluster.cols, name="")

continuous.colours = c("blue", "white", "yellow")
cont.scale = scale_color_gradientn(colours = continuous.colours)
cont.scale2 = scale_fill_gradientn(colours = continuous.colours)

malignancy.cols = setNames(brewer.pal(n=3,"Set2")[1:2], c("Benign", "Malignant"))
location.cols = setNames(c(gg_color_hue(4), alpha("grey50", 0.2)), c("Adrenal", "Extraadrenal", "Head and neck", "Met","Unknown"))

