###########################################################################################
# author: Tanya Karagiannis
# name: seurat_clustering.R
# purpose: script to process and cluster the scRNA-seq of PBMCs (with LPS treatment) of opioid dependent individuals
###########################################################################################

# Load libraries

###########################################################################################

library(dplyr) 
library(Matrix) 
library(vioplot) 
library(useful) 
library(ggplot2)
library(Seurat) #v2.3
library(magrittr)

###########################################################################################

# Analysis

###########################################################################################

#load data from cellranger output
Combined_PBMC <- Read10X("filtered_gene_bc_matrices_mex/GRCh38")



#set seurat object
Combined_PBMC <- CreateSeuratObject(raw.data = Combined_PBMC,
                min.cells=10, min.genes=300,
                project="Combined_PBMC")



# Add mitochondrial expressing genes to metadata

mito.genes <- grep("^mt-", rownames(Combined_PBMC@data), value = T)
percent.mito <- Matrix::colSums(Combined_PBMC@data[mito.genes, ])/Matrix::colSums(Combined_PBMC@data)

Combined_PBMC <- AddMetaData(Combined_PBMC, percent.mito, "percent.mito")
VlnPlot(object = Combined_PBMC, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

par(mfrow = c(1, 2))
GenePlot(object = Combined_PBMC, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = Combined_PBMC, gene1 = "nUMI", gene2 = "nGene")


# Filter gene counts that have unique gene counts over 2000 and less than 300. Filter UMI over 10000. Filter out percent of mitochondrial expressed genes greater than 15%.

Combined_PBMC <- FilterCells(object = Combined_PBMC, subset.names = c("nGene", "nUMI"), low.thresholds = c(300, 0, -Inf), high.thresholds = c(2000,10000, 0.15))

#PCA Analysis

Combined_PBMC <- NormalizeData(object = seurat_object, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableGenes(mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, y.cutoff = 2) %>%
  ScaleData(vars.to.regress= c("nUMI")) %>%
  RunPCA(pc.genes = Combined_PBMC@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50) %>%
  ProjectPCA(object = Combined_PBMC, do.print = FALSE) %>%

#Cluster Analysis
Combined_PBMC <- FindClusters(object = Combined_PBMC_Opioid, reduction.type = "pca", dims.use = 1:20, resolution = 1.5 , print.output = 0, save.SNN = TRUE, force.recalc = TRUE) %>% 
  RunTSNE(dims.use = 1:20)


TSNEPlot(object = Combined_PBMC_Opioid, do.label = TRUE)



#Set cell types

Combined_PBMC_Opioid <- SetAllIdent(object = Combined_PBMC_Opioid, id = "res.1.5")

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, 23,24, 25)
new.cluster.ids <- c("Naive CD4+ T","Naive CD4+ T","Naive CD8+ T","Naive CD4+ T","LPS CD4+ T","Naive CD4+ T", "Naive B Cells","Naive CD8+ T Memory","Naive CD4+ T", "LPS Activated T", "LPS CD8+ T", "LPS CD4+ T", "Naive NK Cells", "Naive B Cells", "LPS B Cells", "Naive Monocytes", "LPS CD4+ T", "LPS CD4+ T","LPS Activated T", "NA","LPS NK Cells", "NA", "LPS Monocytes", "NA","NA","NA") 
Combined_PBMC_Opioid@ident <- plyr::mapvalues(x = Combined_PBMC_Opioid@ident, from = current.cluster.ids, to = new.cluster.ids)

Combined_PBMC_Opioid<- StashIdent(object = Combined_PBMC_Opioid, save.name = "Initial_Cell_Types")



#Remove clusters that have very high UMI counts (labeled NA) 
#Assign major cell types under "Final_Cell_Types" in metadata
PBMC_Final <- SubsetData(Combined_PBMC_Opioid, ident.remove = c("NA"))

PBMC_Final <- SetAllIdent(object = PBMC_Final, id = "res.1.5")

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,22)
new.cluster.ids <- c("CD4+ T Cells","CD4+ T Cells","CD8+ T Cells Naive","CD4+ T Cells","CD4+ T Cells","CD4+ T Cells", "B Cells","CD8+ T Cells Memory","CD4+ T Cells", "CD4+ T Cells", "CD8+ T Cells Memory", "CD4+ T Cells", "NK Cells", "B Cells", "B Cells", "Monocytes", "CD4+ T Cells", "CD4+ T Cells","CD4+ T Cells", "NK Cells", "Monocytes") 
PBMC_Final@ident <- plyr::mapvalues(x = PBMC_Final@ident, from = current.cluster.ids, to = new.cluster.ids)

PBMC_Final<- StashIdent(object = PBMC_Final, save.name = "Final_Cell_Types")


