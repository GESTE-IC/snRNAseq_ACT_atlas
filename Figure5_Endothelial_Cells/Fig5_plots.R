
library(Matrix)
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(sctransform)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(pheatmap)
library(monocle)
library(SeuratWrappers)
library(clusterProfiler)
library(viridis)
library(clustree)

dataset_endoth <- readRDS("/path/to/Seurat_objects/Merge/dataset_endothelial_cells.RDS")
Go_enricher <- readRDS("/path/to/enrichment/results/GO_enricher_endothelial_cellss.RDS")

cds <- readRDS("/path/to/Monocle/results/Cds_Endothelial_Cells.RDS")


### Clustree

clustering_info <- dataset_endoth@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")


### UMAP

my.cell.type.cols <- c("#3E9BFEFF","#30123BFF", "#46F884FF", "#7A0403FF",  "#F05B12FF",  "#E1DD37FF")

p_umap_celltype <- DimPlot(dataset_endoth, reduction = "umap", group.by = "cell_type", label = FALSE, repel = FALSE, raster = FALSE) + 
  scale_color_manual(values = my.cell.type.cols) +
  ggtitle("")


### Stacked Violin Plot for marker genes

features <- c("NRP1", "ADGRL2", "VWF","ANGPT2", "MMRN1","PROX1", "FBLN5","GJA5", "HSPB1", "HSPA1A", "PLAT","KCNIP4")

dataset_endoth$celltype_reorder <- factor(dataset_endoth$new_annot, levels = c("EC-venous", "EC-HSP+","EC-arterial", "EC-lymphatic", "TEC2", "TEC1")) 
 
p_Vln_stacked <- VlnPlot(dataset_endoth, features = features,  cols = my.cell.type.cols, fill.by = 'ident', assay = "RNA", pt.size =, 0, adjust=1.5, stack = TRUE) +
  ggtitle("") +
  scale_y_discrete(limits = levels(dataset_endoth$celltype_reorder)) +
  NoLegend() +
  theme(plot.title = element_text(size = 15, hjust = 0.5), axis.title = element_blank())


### Proportion of cells in each histotype

df <- as.data.frame(table(dataset_endoth$histotype, dataset_endoth$cell_type))
prop_ggplot <- ggplot(df, aes(Var1, Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = c("#3E9BFEFF","#30123BFF", "#46F884FF", "#7A0403FF",  "#F05B12FF",  "#E1DD37FF")) + 
  theme(axis.text.x=element_text(angle=45, size=7, hjust=1),
      strip.background = element_rect(colour="black", fill="white", size=0.5, linetype="solid"),
      strip.text.x = element_text(size=8, face = "bold"),
      axis.text = element_text(size=7), 
      axis.title = element_text(size=8, face = "bold"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.key.size = unit(0.5, 'cm'), 
      legend.key.height = unit(0.5, 'cm'), 
      legend.key.width = unit(0.5, 'cm'), 
      legend.title = element_text(size=7),
      legend.text = element_text(size=7) ) + 
  ylab("Proportion of endothelial cells")  +
  xlab("Sample type") +
  guides(fill=guide_legend(title="Transcriptome cluster"))


### DotPlot Enrichment analyses

dotplot(Go_enricher, showCategory = 4, font.size = 10)

 
### Trajectory DDRTree

p_DDRT_Pseudotime <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = T, show_tree = TRUE, cell_size = 0.5) +  
        scale_color_viridis()

p_DDRT_Celltype <-  plot_cell_trajectory(cds, color_by = "cell_type", cell_size = 0.5) + 
  scale_color_manual(values = c("#46F884FF", "#7A0403FF", "#E1DD37FF",  "#F05B12FF" ,"#3E9BFEFF", "#30123BFF"))


### Trajectory Genes in Pseudotime

my_genes <- c("ANGPT2", "VWF", "KCNQ3", "LAMB1", "ANO2", "ENPP2")

cds_subset <- cds[my_genes,]
p_genes_branch <- plot_genes_branched_pseudotime(cds_subset, branch_point = 2, color_by = "cell_type", ncol = 2, panel_order = c("VWF", "ANGPT2", "ANO2", "KCNQ3", "LAMB1", "ENPP2")) +
  scale_color_manual(values = c( "#46F884FF", "#7A0403FF", "#E1DD37FF",  "#F05B12FF" ,"#3E9BFEFF", "#30123BFF")) 
  
              
