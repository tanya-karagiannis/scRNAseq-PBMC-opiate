###########################################################################################
# author: Tanya Karagiannis
# name: seurat_hashing.R
# purpose: script to process and demultiplex scRNA-seq of IFNb treated PBMCs of opioid-dependent individuals and controls
###########################################################################################

# Load libraries

###########################################################################################

library(dplyr) 
library(Matrix) 
library(vioplot) 
library(useful) 
library(ggplot2)
library(Seurat) #v3

###########################################################################################

# Analysis

###########################################################################################

#load RNA counts

RNA_count <- Read10X(data.dir = "./filtered_feature_bc_matrix/")
PBMC_IFNB <- CreateSeuratObject(counts = RNA_count, project = "PBMC_IFNB", min.cells = 1, min.features = 0)

#load HTO counts
Hashtag_count <- read.csv("./Hashtag_IFNB_Count.txt",header=T,row.names=1)

#intersection between RNA and HTO counts

joint_bcs <- intersect(colnames(PBMC_IFNB[["RNA"]]@data), colnames(Hashtag_count))

length(colnames(PBMC_IFNB[["RNA"]]@data))
length(colnames(Hashtag_count))
length(joint_bcs)


# Subset RNA and HTO counts by joint cell barcodes

PBMC_IFNB <- subset(PBMC_IFNB, cells = c(joint_bcs))

Hashtag_count <- as.matrix(Hashtag_count[,joint_bcs])


###Add HTO counts to the RNA Seurat object metadata.

# Highly variable genes
PBMC_IFNB <- FindVariableFeatures(PBMC_IFNB)
PBMC_IFNB <- ScaleData(PBMC_IFNB)

# Add HTO data as a new assay independent from RNA
PBMC_IFNB[["HTO"]] <- CreateAssayObject(counts = Hashtag_count)



# Normalize HTO data using centered log-ratio (CLR) transformation
PBMC_IFNB <- NormalizeData(PBMC_IFNB, assay = "HTO", normalization.method = "CLR")




###Demultiplex samples based on cell hashing
#Used default settings

PBMC_IFNB <- HTODemux(PBMC_IFNB, assay = "HTO", positive.quantile =  0.99)

#Look at global classification result: singlets, doublets, or negative/ambiguous cells

table(PBMC_IFNB@meta.data$HTO_classification.global)


###Visualize cells based on hto global classification: grouping cells by singlets and doublets
#remove negative cells from the object
PBMC_IFNB.subset <- subset(PBMC_IFNB, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = PBMC_IFNB.subset, assay = "HTO"))))

# tSNE embedding visualization
PBMC_IFNB.subset <- RunTSNE(PBMC_IFNB.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
DimPlot(PBMC_IFNB.subset)

# visualize using heatmap
HTOHeatmap(PBMC_IFNB, assay="HTO",ncells = 5000)

###Calculate percent of mitochondrial genes expressed and filter RNA data
PBMC_IFNB_singlet[["percent.mt"]] <- PercentageFeatureSet(PBMC_IFNB_singlet, pattern = "^MT-")


VlnPlot(PBMC_IFNB_singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

PBMC_IFNB_singlet <- subset(PBMC_IFNB_singlet, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

###Cluster and visualize the final HTO classification of the cells using tSNE.
#Keep only singlets and label combined data based on singlets and sample hto identification

Idents(PBMC_IFNB) <- "HTO_classification.global"

# Extract the singlets
PBMC_IFNB_singlet <- subset(PBMC_IFNB,idents = "Singlet")

Idents(PBMC_IFNB_singlet) <- "HTO_classification"

PBMC_IFNB_singlet$sample_condition <- Idents(PBMC_IFNB_singlet)

all.genes <- rownames(PBMC_IFNB_singlet)

#Normalize Data
PBMC_IFNB_singlet <- NormalizeData(PBMC_IFNB_singlet, normalization.method = "LogNormalize", scale.factor = 10000) %>%
	ScaleData(features = all.genes) %>%
	RunPCA(features = VariableFeatures(PBMC_IFNB))


PBMC_IFNB_singlet <- FindNeighbors(PBMC_IFNB_singlet, reduction = "pca", dims = 1:10) %>%
	FindClusters(resolution = 0.6) %>%
	RunTSNE(reduction = "pca", dims = 1:10)

# Projecting singlet identities on TSNE visualization
DimPlot(PBMC_IFNB_singlet, group.by = "sample_condition", reduction = "tsne")

###Set cell type identities

new.cluster.ids <- c("T Cells","T Cells","T Cells","T Cells", "T Cells", "B Cells","T Cells","NK Cells","T Cells","T Cells","T Cells","Monocytes")
names(new.cluster.ids) <- levels(PBMC_IFNB_singlet)
PBMC_IFNB_singlet <- RenameIdents(PBMC_IFNB_singlet, new.cluster.ids)
DimPlot(PBMC_IFNB_singlet, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

PBMC_IFNB_singlet$Cell_Types <- Idents(PBMC_IFNB_singlet)
