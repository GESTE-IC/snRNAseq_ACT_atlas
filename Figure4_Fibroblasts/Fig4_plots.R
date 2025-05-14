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

dataset_fibro <- readRDS("/path/to/Seurat_objects/Merge/dataset_fibroblasts.RDS")
Go_enricher <- readRDS("/path/to/enrichment/results/GO_enricher_fibroblasts.RDS")

cds <- readRDS("/path/to/Monocle/results/Cds_Fibroblasts.RDS")


### Clustree

clustering_info <- dataset_fibro@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")


### UMAP

my.cell.type.cols <- c("#28BBECFF", "#466BE3FF", "#30123BFF", "#31F299FF", "#A2FC3CFF", "#EDD03AFF", "#FB8022FF", "#D23105FF", "#7A0403FF")

p_umap_celltype <- DimPlot(dataset_fibro, reduction = "umap", group.by = "cell_type", label = FALSE, repel = FALSE, raster = FALSE) + 
  scale_color_manual(values = my.cell.type.cols) +
  ggtitle("")


### Stacked Violin Plot for marker genes

features <- c("FN1", "PDGFRB", "HDAC9", "FGF7", "EFEMP1", "POSTN", "PLA2G5", "GJC1", "ADIPOQ", "PLIN1", "CSPG4", "MYH11", "PRIMA1", "GINS3", "NR4A1", "COL4A4", "COL4A3") 

dataset_fibro$celltype_reorder <- factor(dataset_fibro$cell_type, levels = c("Resident fibroblasts 3", "Resident fibroblasts 2", "Resident fibroblasts 1",
                                                                      "Schwann cells", "Pericytes", "Adipocytes",
 
p_Vln_stacked <- VlnPlot(dataset_fibro, features = features,  cols = my.cell.type.cols, fill.by = 'ident', assay = "RNA", pt.size =, 0, adjust=1.5, stack = TRUE) +
  ggtitle("") +
  scale_y_discrete(limits = levels(dataset_fibro$celltype_reorder)) +
  NoLegend() +
  theme(plot.title = element_text(size = 15, hjust = 0.5), axis.title = element_blank())


### Proportion of cells in each histotype

df <- as.data.frame(table(dataset_fibro$histotype, dataset_fibro$cell_type))
prop_ggplot <- ggplot(df, aes(Var1, Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = c("#28BBECFF", "#466BE3FF", "#30123BFF", "#31F299FF", "#A2FC3CFF", "#EDD03AFF", "#FB8022FF", "#D23105FF", "#7A0403FF")) + 
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
  ylab("Proportion of fibroblasts")  +
  xlab("Sample type") +
  guides(fill=guide_legend(title="Transcriptome cluster"))


### DotPlot Enrichment analyses

dotplot(Go_enricher, showCategory = 4, font.size = 10)

 
### Trajectory DDRTree

p_DDRT_Pseudotime <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = T, show_tree = TRUE, cell_size = 0.5) +  
        scale_color_viridis()

p_DDRT_Celltype <-  plot_cell_trajectory(cds, color_by = "cell_type", cell_size = 0.5) + 
  scale_color_manual(values = c("#28BBECFF", "#466BE3FF", "#30123BFF", "#FB8022FF", "#D23105FF", "#7A0403FF"))


### Trajectory Genes in Pseudotime

my_genes_1 <- c("NGF", "RGS5", "SEMA5A", "CD36")
my_genes_2 <- c("BNC2", "BICC1", "EFEMP1", "VCAN") 

cds_subset_1 <- cds[my_genes_1,]
p_genes_branch_1 <- plot_genes_branched_pseudotime(cds_subset_1, branch_point = 1, color_by = "cell_type", ncol = 4, panel_order = c("RGS5", "SEMA5A", "NGF", "CD36")  ) +
  scale_color_manual(values =  c("#28BBECFF", "#466BE3FF", "#30123BFF", "#FB8022FF", "#D23105FF", "#7A0403FF")) 

cds_subset_2 <- cds[my_genes_2,]
p_genes_branch_2 <- plot_genes_branched_pseudotime(cds_subset_2, branch_point = 1, color_by = "cell_type", ncol = 4, panel_order = c("BNC2", "BICC1", "EFEMP1", "VCAN")  ) +
  scale_color_manual(values =  c("#28BBECFF", "#466BE3FF", "#30123BFF", "#FB8022FF", "#D23105FF", "#7A0403FF")) 
  
              


