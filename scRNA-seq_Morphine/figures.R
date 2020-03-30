###########################################################################################
# author: Tanya Karagiannis
# name: figure_morphine.R
# purpose: script for figures for scRNA-seq of healthy PBMCs treated with morphine in vitro
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

#average geneset expression
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
	
	cell_name <- subset(pbmc_singlet, idents = cell.type[n])
  
  Idents(cell_name) <- "hto_classification"
  
  cell_name <- subset(cell_name, idents = c("Mock+LPS","+Mor+LPS"))
  #cell_name <- subset(cell_name, idents = c("Mock","+Mor","Mock+LPS","+Mor+LPS"))
  
  cell_name@meta.data$hto_classification <- factor(x = cell_name@meta.data$hto_classification, levels = c("Mock+LPS","+Mor+LPS"))
  #cell_name@meta.data$hto_classification <- factor(x = cell_name@meta.data$hto_classification, levels = c("Mock+LPS","+Mor+LPS"))
  

  Idents(cell_name) <- "hto_classification"
  
  cell_name <- ScaleData(cell_name, features = rownames(cell_name))
  
  cell_name <- avgGeneExp(cell_name)

  #my_comparisons <- list(c("Mock+LPS","+Mor+LPS"))
  #my_comparisons <- list(c("Mock", "+Mor","Mock+LPS","+Mor+LPS"))
	
	p <- VlnPlot(object = cell_name, features = geneset, group.by = group,
               x.lab.rot = TRUE,
               cols = c(
                #"#0B6623", "#800080", 
                "#CCCC00", "#1E90FF")
               ) +
    geom_boxplot(width = 0.2) +
    #stat_compare_means(comparisons=my_comparisons, label = "p.signif", method = "t.test") +
    #stat_compare_means(label.y = 2)+
    ylim(0,2)+
    theme(legend.text=element_text(family = "Arial"),
          axis.title.x = element_text(family="Arial"),
          axis.title.y = element_text(family="Arial", size = 30),
          axis.text.x = element_text(family="Arial"),
          axis.text.y = element_text(family="Arial", size = 30)
    )
	}

#Purpose: scaled expression heatmaps for genes in each geneset across all cells
htmap <- function(seurat_object, cell.type = NULL, geneset = Antiviral_Genes, labels = c("ISG15","ISG20","MX1","IFIT1","IFIT2","IFIT3","TNF","NFKB1","TRAF1","CCL5","CCL3","CCL4")){
  	#subset expression data for specific cell type and extract scaled expression matrix
  	cell_name <- subset(seurat_object, idents = cell.type[n])

    cell_name <- ScaleData(cell_name, features = rownames(cell_name))

    Idents(cell_name) <- "hto_classification"

    #Mock <- subset(cell_name, idents = c("Mock"))
    #Mor <- subset(cell_name, idents = c("+Mor"))
    MockLPS <- subset(cell_name, idents = c("Mock+LPS"))
    MorLPS <- subset(cell_name, idents = c("+Mor+LPS"))

    #Data1 <- as.data.frame(as.matrix(Mock[['RNA']]@scale.data))[Antiviral_Genes,]
    #Data2 <- as.data.frame(as.matrix(Mor[['RNA']]@scale.data))[Antiviral_Genes,]
    Data3 <- as.data.frame(as.matrix(MockLPS[['RNA']]@scale.data))[Antiviral_Genes,]
    Data4 <- as.data.frame(as.matrix(MorLPS[['RNA']]@scale.data))[Antiviral_Genes,]


    ha_row = rowAnnotation(df = data.frame(
      Type = factor(c(rep("Core Antiviral",length(Core_Antiviral)),c(rep("Peaked Inflam",length(Peaked_Inflam)), c(rep("Sustained Inflam", length(Sustained_Inflam))))))),
      col = list(Type = c("Core Antiviral" =  "#006400", "Peaked Inflam" = "#DB7093","Sustained Inflam" = "#483D8B")),
      width = unit(10, "mm"),
      annotation_legend_param = list(labels_gp = gpar(fontfamiy= "Arial"))
    )




    #Set samples
    #S1 <- length(Mock@meta.data$hto_classification)
    #S2 <- length(Mor@meta.data$hto_classification)
    S3 <- length(MockLPS@meta.data$hto_classification)
    S4 <- length(MorLPS@meta.data$hto_classification)


    #Heatmaps
    ht1 <- Heatmap(as.matrix(Data1),
                   name = "Scaled Expression",
                   col = colorRamp2(c(-2,0,2), c("purple","black", "yellow")),
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   #top_annotation = ha_column1,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   row_split = c(rep("Core Antiviral",length(Core_Antiviral)),rep("Peaked Inflam",length(Peaked_Inflam)), rep("Sustained Inflam", length(Sustained_Inflam))),
                   column_split = c(rep("Untreated", S1)),
                   column_title_gp = gpar(fill = c("#0B6623"), font = 2)
    )
    ht2 <- Heatmap(as.matrix(Data2),
                   name = "Scaled Expression",
                   col = colorRamp2(c(-2,0,2), c("purple","black", "yellow")),
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   #top_annotation = ha_column1,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   row_split = c(rep("Core Antiviral",length(Core_Antiviral)),rep("Peaked Inflam",length(Peaked_Inflam)), rep("Sustained Inflam", length(Sustained_Inflam))),
                   column_split = c(rep("Morphine", S2)),
                   column_title_gp = gpar(fill = c("#800080"), font = 2)
    )
    # ht3 <- Heatmap(as.matrix(Data3),
    #                name = "Scaled Expression",
    #                col = colorRamp2(c(-2,0,2), c("purple","black", "yellow")),
    #                cluster_rows = FALSE,
    #                cluster_columns = FALSE,
    #                #top_annotation = ha_column1,
    #                show_row_names = FALSE,
    #                show_column_names = FALSE,
    #                row_split = c(rep("Core Antiviral",length(Core_Antiviral)),rep("Peaked Inflam",length(Peaked_Inflam)), rep("Sustained Inflam", length(Sustained_Inflam))),
    #                column_split = c(rep("LPS", S3)),
    #                column_title_gp = gpar(fill = c("#CCCC00"), font = 2)
    # )
    # ht4 <- Heatmap(as.matrix(Data4),
    #                name = "Scaled Expression",
    #                col = colorRamp2(c(-2,0,2), c("purple","black", "yellow")),
    #                cluster_rows = FALSE,
    #                cluster_columns = FALSE,
    #                #top_annotation = ha_column1,
    #                show_row_names = FALSE,
    #                show_column_names = FALSE,
    #                row_split = c(rep("Core Antiviral",length(Core_Antiviral)),rep("Peaked Inflam",length(Peaked_Inflam)), rep("Sustained Inflam", length(Sustained_Inflam))),
    #                column_split = c(rep("Morphine + LPS", S4)),
    #                column_title_gp = gpar(fill = c("#1E90FF"), font = 2)
    # )

    #cell <- gsub(" ","_",cell_types[n])

    loc.labels <- which(Antiviral_Genes %in% labels)
    labels <- Antiviral_Genes[loc.labels]

    #pdf(paste0("./PBMC_", cell, "_Antiviral.pdf"), 10,5)
    draw(ha_row + ht3 + ht4 + rowAnnotation(link = anno_mark(at = loc.labels, labels = labels, which = c("row"))))
    #dev.off()

}

###########################################################################################
#Figures
###########################################################################################

#Figure 3c

cell_types <- c("CD4 T Cells", "CD8 T Cells","B Cells", "NK Cells", "Monocytes") 

#Heatmap of scaled expression for genesets across sample cells in each cell type
lapply(1:length(cell_types), function(x){
	p <- htmap(PBMC_Morphine, cell.type = cell_types[x])
	})

#Violin plots of average geneset expression across sample cells in each cell type
lapply(1:length(cell_types), function(x){
	p <- vlnplot(PBMC_Morphine, cell.type = cell_types[x], geneset = c("core","peaked","sustained"), group = "sample_condition")
	})

#Supplementary Figure 16


#tSNE embedding with cell type labels across all cells
DimPlot(PBMC_Morphine, group.by = "Cell_Types",reduction = "tsne")

#tSNE embedding with sample labels across all cells
sample_colors <- c("#0B6623", "#800080","#CCCC00", "#1E90FF")
DimPlot(PBMC_Morphine, group.by = "sample_condition", cols = sample_colors, reduction = "tsne")

#expression of canonical gene markers for each cell type in tSNE embedding
FeaturePlot(object = PBMC_Morphine, features = c("CD3D","IL7R","CD8A","CD8B","MS4A1","CD79","GNLY","NKG7","LYZ","CD14","FCER1A","CST3"), cols = c("grey", "darkred"))


#Supplement Figures 17-21
#uncomment sections in function htmap and vlnplot to include naive samples (Mock and +Mor) and LPS samples (Mock+LPS and +Mor+LPS)

cell_types <- c("CD4 T Cells", "CD8 T Cells","B Cells", "NK Cells", "Monocytes") 

#Heatmap of scaled expression for genesets across sample cells in each cell type
lapply(1:length(cell_types), function(x){
  p <- htmap(PBMC_Morphine, cell.type = cell_types[x])
  })

#Violin plots of average geneset expression across sample cells in each cell type
lapply(1:length(cell_types), function(x){
  p <- vlnplot(PBMC_Morphine, cell.type = cell_types[x], geneset = c("core","peaked","sustained"), group = "sample_condition")
  })

