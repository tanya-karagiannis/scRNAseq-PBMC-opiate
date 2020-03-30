###########################################################################################
# author: Tanya Karagiannis
# name: seurat_hashing.R
# purpose: script to process and demultiplex scRNA-seq of healthy PBMCs treated with morphine in vitro
###########################################################################################

# Load libraries

###########################################################################################

library(dplyr) 
library(Matrix) 
library(vioplot) 
library(useful) 
library(ggplot2)
library(Seurat) #v2.3

###########################################################################################

# Analysis

###########################################################################################

#load RNA counts
PBMC_Morphine <- Read10X("filtered_gene_bc_matrices/GRCh38/")

PBMC_Morphine <- CreateSeuratObject(raw.data = PBMC_Morphine,
                min.cells=10, min.genes=0,
                project="Combined PBMC Morphine")

#load HTO counts
PBMC_HTO=read.csv("PBMC_HTO.txt",header=T,row.names=1)


#intersection between RNA and HTO counts

joint_bcs <- intersect(colnames(PBMC_Morphine@data),colnames(PBMC_HTO))

# Subset RNA and HTO counts by joint cell barcodes

PBMC_Morphine <- SubsetData(PBMC_Morphine, cells.use = joint_bcs)

PBMC_HTO <- as.matrix(PBMC_HTO[,joint_bcs])

# Confirm that the HTO have the correct names
print (colnames(PBMC_HTO))



###Add HTO counts to the RNA Seurat object metadata

# Highly variable genes
PBMC_Morphine <- FindVariableGenes(PBMC_Morphine,do.plot = F,display.progress = FALSE)
PBMC_Morphine <- ScaleData(PBMC_Morphine, genes.use = PBMC_Morphine@var.genes,display.progress = FALSE, num.cores=4)

# Add HTO data as a new assay independent from RNA
PBMC_Morphine <- SetAssayData(PBMC_Morphine,assay.type = "HTO",slot = "raw.data",new.data = PBMC_HTO)



### Normalize HTO data using centered log-ratio (CLR) transformation
PBMC_Morphine <- NormalizeData(PBMC_Morphine,assay.type = "HTO",normalization.method = "genesCLR",display.progress = FALSE)


### Demultiplex samples based on cell hashing
#Using default settings

PBMC_Morphine <- HTODemux(PBMC_Morphine,assay.type = "HTO",positive_quantile =  0.99,print.output = FALSE, k_function = "clara")


#Look at global classification result: shows assigned singlets, doublets, or negative/ambiguous cells
print(table(PBMC_Morphine@meta.data$hto_classification_global))



###Visualize cells based on HTO global classification: grouping cells by singlets and doublets.
#Supplementary Figures


#Remove negative cells from the object
PBMC_Morphine <- SetAllIdent(PBMC_Morphine,"hto_classification_global")
PBMC_Morphine_subset <- SubsetData(PBMC_Morphine,ident.remove = "Negative")

# Calculate a distance matrix using HTO
hto_dist_mtx <- as.matrix(dist(t(GetAssayData(PBMC_Morphine_subset,assay.type = "HTO",slot = "data"))))

# tSNE embedding visualization
PBMC_Morphine_subset <- RunTSNE(PBMC_Morphine_subset, distance.matrix = hto_dist_mtx,perplexity=100)
TSNEPlot(PBMC_Morphine_subset,pt.size = 0.2)

# visualize using heatmap
HTOHeatmap(PBMC_Morphine,num.cells = 5000)


###Cluster and visualize the final HTO classification of the cells using tSNE.
#Keep only singlets for this process

# Extract the singlets
PBMC_Morphine <- SetAllIdent(PBMC_Morphine,id = "hto_classification_global")

###Cluster Analysis
pbmc_singlet <- SubsetData(PBMC_Morphine,ident.use = "Singlet") %>%
pbmc_singlet <- FindVariableGenes(pbmc_singlet,do.plot = F) %>%
	ScaleData(genes.use = pbmc_singlet@var.genes,display.progress = FALSE) %>%
	RunPCA(pbmc_singlet,pc.genes = pbmc_singlet@var.genes,pcs.print = 0)

pbmc_singlet <- FindClusters(pbmc_singlet,reduction.type = "pca",dims.use = 1:10,print.output = F,resolution = 0.6) %>%
	RunTSNE(reduction.use = "pca",dims.use = 1:10)

# Projecting singlet identities on TSNE visualization
TSNEPlot(pbmc_singlet,group.by = "hto_classification")


#set cell types

pbmc_singlet <- SetAllIdent(object = pbmc_singlet, id = "res.0.6")

current.cluster.ids <- c(0,1,2,3,4,5,6,7,8,9,10)
new.cluster.ids <- c("CD4 T Cells","CD4 T Cells","CD8 T Cells","CD4 T Cells","CD4 T Cells","Monocytes", "Monocytes","B Cells","CD8 T Cells", "NK Cells","NA") 
pbmc_singlet@ident <- plyr::mapvalues(x = pbmc_singlet@ident, from = current.cluster.ids, to = new.cluster.ids)

pbmc_singlet<- StashIdent(object = pbmc_singlet, save.name = "Cell_Types")



