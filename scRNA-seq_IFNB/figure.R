###########################################################################################
# author: Tanya Karagiannis
# name: figure_morphine.R
# purpose: script for figures for scRNA-seq of IFNb-treated PBMCs of opioid-dependent individuals and controls
#generate average expres
###########################################################################################
# Load libraries
###########################################################################################

library(dplyr) 
library(Matrix) 
library(useful) 
library(ggplot2)
library(Seurat) #v3
library(magrittr)


###########################################################################################
#Functions
###########################################################################################

#average expression for all genes in a geneset for each single cell in a cell type 
#Genesets: Three LPS stimulated antiviral genesets: core antiviral, peaked inflam, and sustained inflam
avgGeneExp <- function(seurat_object){
  mean.core <- colMeans(x = as.matrix(seurat_object[['RNA']]@data[Core_Antiviral, ]), na.rm = TRUE)
  mean.peaked <- colMeans(x = as.matrix(seurat_object[['RNA']]@data[Peaked_Inflam, ]), na.rm = TRUE)
  mean.sustained <- colMeans(x = as.matrix(seurat_object[['RNA']]@data[Sustained_Inflam, ]), na.rm = TRUE)
  
  # Add mean expression values
  seurat_object@meta.data$core <- mean.core
  seurat_object@meta.data$peaked <- mean.peaked
  seurat_object@meta.data$sustained <- mean.sustained
}


#violinplot of average geneset expression
vlnplot <- function(seurat_object, cell.type = NULL ,geneset = c("core","peaked","sustained"), group = "sample_condition"){
  
  cell_name <- subset(PBMC_IFNB_singlet, idents = cell.type[n])
  
  cell_name@meta.data$sample_condition <- factor(cell_name@meta.data$sample_condition, 
                                                levels = c("IFNb-C1","IFNb-C2","IFNb-C3","IFNb-O1","IFNb-O2","IFNb-O3"))
  
  Idents(cell_name) <- "sample_condition"
  
  cell_name <- ScaleData(cell_name, features = rownames(cell_name))
  
  
  cell_name <- avgGeneExp(cell_name)

  #my_comparisons <- list(c("Control","Opioid"))
  
  p <- VlnPlot(object = cell_name, features = geneset, group.by = group,
               x.lab.rot = TRUE,
               cols = c("#98FB98", "#4F7942", "#0B6623", "#EE82EE", "#9932CC", "#800080")
  ) +
    geom_boxplot(width = 0.2) +
    #stat_compare_means(comparisons=my_comparisons, label = "p.signif", method = "t.test") +
    #stat_compare_means(label.y = 5)+
    ylim(0,1.5)+
    theme(legend.text=element_text(family = "Arial"),
          axis.title.x = element_text(family="Arial"),
          axis.title.y = element_text(family="Arial", size = 30),
          axis.text.x = element_text(family="Arial"),
          axis.text.y = element_text(family="Arial", size = 30)
    )
  }


###########################################################################################
#Figures
###########################################################################################

#Supplement Figure 14

#tSNE embedding with cell type labels across all cells

DimPlot(PBMC_IFNB_singlet, group.by = "Cell_Types",reduction = "tsne")

#tSNE embedding with sample labels across all cells
Control_Colors <- c("#98FB98", "#4F7942", "#0B6623")
Dependent_Colors <- c("#EE82EE", "#9932CC", "#800080")
DimPlot(PBMC_IFNB_singlet, group.by = "sample_condition", cols = c(Control_Colors, Dependent_Colors),reduction = "tsne")

#expression of canonical gene markers for each cell type in tSNE embedding
FeaturePlot(object = PBMC_IFNB_singlet, features = c("CD3D","IL7R","CD8A","CD8B","MS4A1","CD79","GNLY","NKG7","LYZ","CD14","FCER1A","CST3"), cols = c("grey", "darkred"))

#Supplement Figure 15

cell_types <- c("T Cells", "B Cells","NK Cells","Monocytes") 

#Violin plots of average geneset expression across sample cells in each cell type
lapply(1:length(cell_types), function(x){
  p <- vlnplot(PBMC_Morphine, cell.type = cell_types[x], geneset = c("core","peaked","sustained"), group = "sample_condition")
  })
