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

dataset_macro <- readRDS("/path/to/Seurat_objects/Merge/dataset_myeloid_cells.RDS")
Go_enricher <- readRDS("/path/to/enrichment/results/GO_enricher_myeloid_cellss.RDS")

cds <- readRDS("/path/to/Monocle/results/Cds_Myeloid_Cells.RDS")


### Clustree

clustering_info <- dataset_macro@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")


### UMAP

my.cell.type.cols <-  c("#30123BFF", "#4686FBFF","#1AE4B6FF", "#A2FC3CFF",  "#FABA39FF", "#E4460AFF", "#7A0403FF" )

p_umap_celltype <- DimPlot(dataset_macro, reduction = "umap", group.by = "cell_type", label = FALSE, repel = FALSE, raster = FALSE) + 
  scale_color_manual(values = my.cell.type.cols) +
  ggtitle("")


### Stacked Violin Plot for marker genes

features<- c("MRC1", "CD163", "FKBP5", "HS3ST2", "PPARG", "SPP1", "RGS1", "SELENOP","LYVE1", "C3","CX3CR1", "MKI67","TOP2A")

dataset_macro$celltype_reorder <- factor(dataset_macro$cell_type, levels = c("TAM1", "TAM2", "Perivascular TAM", "Inflammatory macrophages",
                                      "Resident macrophages 1", "Resident macrophages 2", "Cycling macrophages")) 
 
p_Vln_stacked <- VlnPlot(dataset_macro, features = features,  cols = my.cell.type.cols, fill.by = 'ident', assay = "RNA", pt.size =, 0, adjust=1, stack = TRUE) +
  ggtitle("") +
  scale_y_discrete(limits = levels(dataset_macro$celltype_reorder)) +
  NoLegend() +
  theme(plot.title = element_text(size = 15, hjust = 0.5), axis.title = element_blank())


### Proportion of cells in each histotype

df <- as.data.frame(table(dataset_macro$histotype, dataset_macro$celltype_reorder))
prop_ggplot <- ggplot(df, aes(Var1, Freq, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_fill_manual(values = my.cell.type.cols) + 
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
  ylab("Proportion of myeloid cells")  +
  xlab("Sample type") +
  guides(fill=guide_legend(title="Transcriptome cluster"))


### DotPlot Enrichment analyses

dotplot(Go_enricher, showCategory = 3, font.size = 10)

 
### Trajectory DDRTree

p_DDRT_Pseudotime <- plot_cell_trajectory(cds, color_by = "Pseudotime", show_branch_points = T, show_tree = TRUE, cell_size = 0.5) +  
        scale_color_viridis()

p_DDRT_Celltype <-  plot_cell_trajectory(cds, color_by = "cell_type", cell_size = 0.5) + 
  scale_color_manual(values = c("#FABA39FF", "#E4460AFF", "#7A0403FF" , "#30123BFF", "#4686FBFF","#1AE4B6FF", "#A2FC3CFF"))


### Trajectory Genes in Pseudotime

my_genes_branch1 <- c("CCL18", "PPARG", "MMP19", "GPNMB", "ABCA1",  "ABCG1") 
cds_subset_1 <- cds[my_genes_branch1,]
p_genes_branch1 <- plot_genes_branched_pseudotime(cds_subset_1, branch_point = 2, color_by = "cell_type", ncol = 2, panel_order = c("PPARG","GPNMB", "CCL18", "MMP19", "ABCA1", "ABCG1")) +
  scale_color_manual(values = c("#FABA39FF", "#E4460AFF", "#7A0403FF" , "#30123BFF", "#4686FBFF","#1AE4B6FF", "#A2FC3CFF")) 

my_genes_branch2 <- c("ADGRB3", "PCNX2", "SORL1", "CX3CR1", "C3",  "PRKG1") 
cds_subset_2 <- cds[my_genes_branch2,]
p_genes_branch2 <- plot_genes_branched_pseudotime(cds_subset_2, branch_point = 2, color_by = "cell_type", ncol = 2, panel_order = c("ADGRB3","PRKG1", "PCNX2", "SORL1", "CX3CR1", "C3")) +
  scale_color_manual(values = c("#FABA39FF", "#E4460AFF", "#7A0403FF" , "#30123BFF", "#4686FBFF","#1AE4B6FF", "#A2FC3CFF)) 
              
