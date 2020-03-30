###########################################################################################
# author: Tanya Karagiannis
# name: seurat_clustering.R
# purpose: script for figures for scRNA-seq of PBMCs (with LPS treatment) of opioid dependent individuals and controls
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

#differential gene expression plots
#differential expresions results from Source Data for each cell type
#interferon genes to label

DE.plot <- function(MAST_res = DE_results, interferon_genes = NULL){

MAST_res$Significance <- ifelse(MAST_res$p_val_adj < 0.05 & abs(MAST_res$avg_logFC) > log2(1.5) & MAST_res$gene %in% interferon_genes, "significant antiviral genes", ifelse(MAST_res$p_val_adj < 0.05 & abs(MAST_res$avg_logFC) > log2(1.5), "significant", "not significant"))

p <- ggplot(MAST_res, aes(x = avg_logFC, y = -log10(p_val_adj))) +
  xlab("log2 fold change") +
  ylab("-log10(FDR)")+
  xlim(-2,2) +
  geom_point(aes(color = Significance), size = 6) +
  scale_color_manual(values = c("grey", "black", "green")) +
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)), linetype="dashed", color = "black", size=2) +
  geom_hline(yintercept = c(-log10(0.05)), linetype="dashed", color = "black", size=2) +
  geom_text_repel(
    data = subset(MAST_res, p_val_adj < 0.05 & abs(avg_logFC) > log2(1.5)),
    aes(label = gene),
    size = 18,
    family = "Arial",
    box.padding = 0.25,
    point.padding = 0.3,
    position = "dodge"
  ) +
  theme_bw() +
  theme(legend.title=element_text(family = "Arial", face = "bold", size=50),
        legend.text=element_text(family = "Arial", size = 50), 
        axis.title.x = element_text(family="Arial", size = 50), 
        axis.title.y = element_text(family="Arial", size = 50),
        axis.text.x = element_text(family="Arial", size = 50),
        axis.text.y = element_text(family="Arial", size = 50), 
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")
  )
}

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
	
	cell_name <- subset(PBMC_Final, idents = cell_types[n])
  
  	Idents(cell_name) <- "ConditionStim"
  	cell_name <- subset(cell_name, idents = c("Control_LPS","Dependent_LPS"))
  
  	cell_name <- ScaleData(cell_name, features = rownames(cell_name))

  	Idents(cell_name) <- "sample_condition"
	
	p <- VlnPlot(object = cell_name, features = geneset, group.by = group,
               x.lab.rot = TRUE,
               cols = c("#EEDD82", "#FFFF33", "#CCCC00",
                            "#87CEFA", "#00BFFF","#1E90FF")
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
  	cell_name <- subset(PBMC_Final, idents = cell_types[n])

	cell_name <- ScaleData(cell_name, features = rownames(cell_name))

	Idents(cell_name) <- "ConditionStim"

	Control <- subset(cell_name, idents = c("Control_LPS"))
	Dependence <- subset(cell_name, idents = c("Dependent_LPS"))

	Data1 <- as.data.frame(as.matrix(Control[['RNA']]@scale.data))[Antiviral_Genes,]
	Data2 <- as.data.frame(as.matrix(Dependence[['RNA']]@scale.data))[Antiviral_Genes,]

	#row annotation
  	ha_row = rowAnnotation(df = data.frame(
  	Type = factor(c(rep("Core Antiviral",length(Core_Antiviral)),c(rep("Peaked Inflam",length(Peaked_Inflam)), c(rep("Sustained Inflam", length(Sustained_Inflam))))))),
 	col = list(Type = c("Core Antiviral" =  "#006400", "Peaked Inflam" = "#DB7093","Sustained Inflam" = "#483D8B")),
  	width = unit(10, "mm"),
  	annotation_legend_param = list(labels_gp = gpar(fontfamiy= "Arial"))
  	)

  	#column annotation
	C1 <- length(Control@meta.data$sample_condition[Control@meta.data$sample_condition == "LPS_Control_3"])
	C2 <- length(Control@meta.data$sample_condition[Control@meta.data$sample_condition == "LPS_Control_4"])
	C3 <- length(Control@meta.data$sample_condition[Control@meta.data$sample_condition == "LPS_Control_5"])
	O1 <- length(Dependence@meta.data$sample_condition[Dependence@meta.data$sample_condition == "LPS_Dependent_10"])
	O2 <- length(Dependence@meta.data$sample_condition[Dependence@meta.data$sample_condition == "LPS_Dependent_11"])
	O3 <- length(Dependence@meta.data$sample_condition[Dependence@meta.data$sample_condition == "LPS_Dependent_12"])
`


	ht1 <- Heatmap(as.matrix(Data1),
               name = "Expression",
               col = colorRamp2(c(-2,0,2), c("purple","black", "yellow")),
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               #top_annotation = ha_column1,
               show_row_names = FALSE,
               show_column_names = FALSE,
               row_split = c(rep("Core Antiviral",length(Core_Antiviral)),rep("Peaked Inflam",length(Peaked_Inflam)), rep("Sustained Inflam", length(Sustained_Inflam))),
               column_split = c(rep("C1", C1), rep("C2", C2), rep("C3", C3)),
               column_title_gp = gpar(fill = c("#EEDD82", "#FFFF33", "#CCCC00"), font = 2)
               )
	ht2 <-  Heatmap(as.matrix(Data2),
                name = "Expression",
                col = colorRamp2(c(-2,0,2), c("purple","black", "yellow")),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                #top_annotation = ha_column2,
                show_row_names = FALSE,
                show_column_names = FALSE,
                row_split = c(rep("Core Antiviral",length(Core_Antiviral)),rep("Peaked Inflam",length(Peaked_Inflam)), rep("Sustained Inflam", length(Sustained_Inflam))),
                column_split = c(rep("O1", O1), rep("O2", O2), rep("O3", O3)),
               column_title_gp = gpar(fill = c("#87CEFA", "#00BFFF", "#1E90FF"), font = 2)
               )

	#cell <- gsub(" ","_",cell_types[n])


	loc.labels <- which(Antiviral_Genes %in% labels)
	labels <- Antiviral_Genes[loc.labels]

	#pdf(paste0("./PBMC_", cell, "_Antiviral.pdf"), 20,10)
	draw(ha_row + ht1 + ht2 + rowAnnotation(link = anno_mark(at = loc.labels, labels = labels, which = c("row"))))
	#dev.off()
}

###########################################################################################
#Figures
###########################################################################################

#Figure 1b

#tSNE of cell types
TSNEPlot(object = PBMC_LPS, do.label = TRUE, pt.size = 0.1, group.by = "Cell_Type", label.size = 10, do.return = TRUE,
  colors.use = colorspace::rainbow_hcl(12))+
  scale_colour_manual(values = colorspace::rainbow_hcl(12), breaks = c(1,2,3,4,5,6,7,8,9,10,11,12), 
    labels = c("(1) Naive CD4+ T","(2) Naive CD8+ T","(3) Naive CD8+ T Memory","(4) Naive B Cells","(5) Naive NK Cells","(6) Naive Monocytes","(7) LPS CD4+ T", "(8) LPS CD8+ T","(9) LPS Activated T","(10) LPS B Cells","(11) LPS NK Cells", "(12) LPS Monocytes")) +
  guides(colour = guide_legend(override.aes = list(alpha = 0))) +
  theme(legend.text=element_text(family = "Arial", size = 30), 
    axis.title.x = element_text(family="Arial", size = 30), 
    axis.title.y = element_text(family="Arial", size = 30),
    axis.text.x = element_text(family="Arial", size = 30),
    axis.text.y = element_text(family="Arial", size = 30)
    )



Control_Colors <- c("#D0F0C0","#98FB98","#50C878","#3BB143","#00A86B","#4F7942","#0B6623")
Dependent_Colors <- c("#D8BFD8","#EE82EE", "#FF00FF", "#BA55D3","#8A2BE2", "#9932CC", "#800080")
LPS_Colors <- c("#EEDD82","#FFFF33", "#CCCC00", "#87CEFA", "#00BFFF","#1E90FF")

#tSNE of samples
TSNEPlot(object = PBMC_Final, do.label = FALSE, pt.size = 0.1, group.by = "sample_condition", 
  plot.order = rev(c(paste("Control_", 1:7,sep=""),paste("Dependent_", 8:14,sep=""), paste("LPS_Control_", 3:5,sep=""),paste("LPS_Dependent_", 10:12,sep=""))),
  colors.use = c(Control_Colors, Dependent_Colors, LPS_Colors),
  do.return = TRUE) + 
  theme(legend.text=element_text(family = "Arial", size = 25), 
    axis.title.x = element_text(family="Arial", size = 30), 
    axis.title.y = element_text(family="Arial", size = 30),
    axis.text.x = element_text(family="Arial", size = 30),
    axis.text.y = element_text(family="Arial", size = 30)
    )

#Figure 1c

DE.plot(MAST_Mono_Naive, interferon_genes = c("IFIT2","PMAIP1","ISG15","IFIT3","ZC3HAV1","GBP2","GBP5"))
DE.plot(MAST_CD4_LPS, interferon_genes = c("ISG15","IFIT3","IFIT2","OASL","IFIT1","MX1","PMAIP1","LY6E","RSAD2","ISG20","IFI6","IFI44"))
DE.plot(MAST_T_Activated, interferon_genes = c("IFIT1","IFIT2","IFIT3","OASL","PMAIP1","ISG15","IFI44","MX1"))
DE.plot(MAST_NK_LPS, interferon_genes = c(c("IFIT2","PMAIP1","ISG15","IFIT3","ZC3HAV1","GBP2","GBP5")))


#Figure 2

cell_types <- c("LPS CD4+ T", "LPS Activated T","LPS B Cells","LPS CD8+ T","LPS NK Cells", "LPS Monocytes") 

#Heatmap of scaled expression for genesets across control cells and opioid cells in each cell type
lapply(1:length(cell_types), function(x){
	p <- htmap(PBMC_LPS, cell.type = cell_types[x])
	})

#Violin plots of average geneset expression across control and opioid cell in each cell type
lapply(1:length(cell_types), function(x){
	p <- vlnplot(PBMC_LPS, cell.type = cell_types[x], geneset = c("core","peaked","sustained"), group = "sample_condition")
	})
