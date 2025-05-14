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

dataset_lympho <- readRDS("/path/to/Seurat_objects/Merge/dataset_lymphoid_cells.RDS")
Go_enricher <- readRDS("/path/to/enrichment/results/GO_enricher_lymphoid_cellss.RDS")

cds <- readRDS("/path/to/Monocle/results/Cds_Lymphoid_Cells.RDS")


### Clustree

clustering_info <- dataset_lympho@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")


### UMAP

my.cell.type.cols <-  c("#7A0403FF", "#FB8022FF", "#A2FC3CFF",  "#28BBECFF", "#30123BFF")

p_umap_celltype <- DimPlot(dataset_lympho, reduction = "umap", group.by = "cell_type", label = FALSE, repel = FALSE, raster = FALSE) + 
  scale_color_manual(values = my.cell.type.cols) +
  ggtitle("")


### Stacked Violin Plot for marker genes

features <- c("IGLL5", "FCRL5", "IL7R", "RUNX2", "GNLY", "KLRF1", "TOX", "CCL5")

dataset_lympho$celltype_reorder <- factor(dataset_lympho$new_annot, levels = c("Exhausted T cells","Unassigned T cells", "NK-like T cells", "Naive/Memory T cells", "Plasma B cells")) 
 
p_Vln_stacked <- VlnPlot(dataset_lympho, features = features,  cols = rev(my.cell.type.cols), fill.by = 'ident', assay = "RNA", pt.size =, 0, adjust=1, stack = TRUE) +
  ggtitle("") +
  scale_y_discrete(limits = levels(dataset_lympho$celltype_reorder)) +
  NoLegend() +
  theme(plot.title = element_text(size = 15, hjust = 0.5), axis.title = element_blank())


### Proportion of cells in each histotype

df <- as.data.frame(table(dataset_lympho$histotype, dataset_lympho$cell_type))
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
  ylab("Proportion of lymphoid cells")  +
  xlab("Sample type") +
  guides(fill=guide_legend(title="Transcriptome cluster"))


### DotPlot Enrichment analyses

dotplot(Go_enricher, showCategory = 3, font.size = 10)

 
### Trajectory DDRTree

p_DDRT_Pseudotime <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = T, show_tree = TRUE, cell_size = 0.5) +  
        scale_color_viridis()

p_DDRT_Celltype <-  plot_cell_trajectory(cds, color_by = "cell_type", cell_size = 0.5) + 
  scale_color_manual(values = c("#30123BFF",  "#28BBECFF", "#A2FC3CFF", "#FB8022FF"))


### Trajectory Genes in Pseudotime

my_genes <- c("THEMIS", "INPP4B", "KLRF1", "KLRD1", "GNLY",  "BNC2") 

cds_subset <- cds[my_genes,]
p_genes_branch <- plot_genes_branched_pseudotime(cds_subset, color_by = "cell_type", ncol = 2, panel_order = c("KLRF1", "KLRD1", "GNLY", "BNC2", "INPP4B", "THEMIS")) +
  scale_color_manual(values = c( "#46F884FF", "#7A0403FF", "#E1DD37FF",  "#F05B12FF" ,"#3E9BFEFF", "#30123BFF")) 
  
              
