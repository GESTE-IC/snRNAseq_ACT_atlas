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

features <- c("FN1", "PDGFRB", "HDAC9", "FGF7", "EFEMP1", "POSTN", "PLA2G5", "GJC1", "ADIPOQ","PLIN1", "CSPG4","MYH11", "PRIMA1","GINS3", "NR4A1", "COL4A4","COL4A3") 

dataset_fibro$celltype_reorder <- factor(dataset_fibro$new_annot, levels = c("Resident fibroblasts 3", "Resident fibroblasts 2", "Resident fibroblasts 1",
                                                                      "Schwann cells", "Pericytes", "Adipocytes",
 
p_Vln_stacked <- VlnPlot(dataset_fibro, features = features,  cols = my.cell.type.cols, fill.by = 'ident', assay = "RNA", pt.size =, 0, adjust=1.5, stack = TRUE) +
  ggtitle("") +
  scale_y_discrete(limits = levels(dataset$celltype_reorder)) +
  NoLegend() +
  theme(plot.title = element_text(size = 15, hjust = 0.5), axis.title = element_blank())


### Proportion of cells in each histotype

df<- as.data.frame(table(dataset_fibro$histotype, dataset_fibro$cell_type))
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


                                                                             
### Trajectory DDRTree -------------------------------------------------


library(Matrix)
library(dplyr)
library(Seurat)
#library(SeuratData)
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

dataset.combined <- readRDS("/path/to/Seurat_objects/Merge/dataset_normal_adrenal_integrated.RDS")
dataset_normaladrenal <- readRDS("/path/to/Seurat_objects/Merge/dataset_normal_adrenal_not_integrated.RDS")
reactome_enricher <- readRDS("/path/to/enrichment/results/Reactome_enricher_normal_adrenal.RDS")

cds <- readRDS("/path/to/Monocle/results/Cds_NormalAdrenal.RDS")
diff_test_res_all <- readRDS("/path/to/Monocle/results/DEG_byPseudotime_NormalAdrenal.RDS")


### Clustree

clustering_info <- dataset.combined@meta.data
clustree(clustering_info, prefix="integrated_snn_res.")


### UMAP

my.cell.type.cols <- c("ZG cells - CYP11B2+" = "#7A0403FF", "ZG cells" = "#F05B12FF", "Intermediate-state cells" = "#E1DD37FF", 
                       "ZF cells" = "#46F884FF",  "ZR cells" = "#3E9BFEFF", "HSP+ cells" = "#30123BFF")

dataset.combined$celltype_integrated <- dataset.combined@active.ident
levels(dataset.combined) <- c("ZG cells - CYP11B2+", "ZG cells", "Intermediate-state cells", "ZF cells", "ZR cells", "HSP+ cells")
dataset_normaladrenal$celltype_integrated <- dataset.combined@active.ident

p_umap_celltype <- DimPlot(dataset.combined, reduction = "umap", group.by = "celltype_integrated", label = FALSE, raster = FALSE, cols = my.cell.type.cols)

p_umap_celltype_not_integrated <- DimPlot(dataset_normaladrenal, reduction = "umap", group.by = "celltype_integrated", label = FALSE, raster = FALSE, cols = my.cell.type.cols)

p_umap_samples_not_integrated <- DimPlot(dataset_normaladrenal, reduction = "umap", group.by = "orig.ident", label = FALSE, raster = FALSE) + 
  scale_color_viridis_d(option = "plasma", direction = -1)


### Stacked Violin Plot for marker genes

features <- c("CYP11B2", "DACH1", "ANO4", "HSD3B2","CYP11B1","CYP17A1", "SULT2A1","HSP90AA1")

dataset.combined$celltype_reorder <- factor(dataset.combined$celltype_integrated, 
                                            levels = c("HSP+ cells", "ZR cells", "ZF cells", "Intermediate-state cells", "ZG cells", "ZG cells - CYP11B2+"))

p_Vln_stacked <- VlnPlot(dataset.combined, features = features,  cols = my.cell.type.cols, fill.by = 'ident', assay = "RNA", pt.size =, 0, adjust=1.5, stack = TRUE) +
  ggtitle("") +
  scale_y_discrete(limits = levels(dataset$celltype_reorder)) +
  NoLegend() +
  theme(plot.title = element_text(size = 15, hjust = 0.5), axis.title = element_blank())


### Proportion of cells in each sample

temp_Names <- unique(dataset.combined$orig.ident)
dataset@meta.data$cell_type <- dataset.combined@active.ident
temp_colname <- "cell_type"
colour_pal <- data.frame(celltype = c("ZG cells - CYP11B2+","ZG cells","Intermediate-state cells","ZF cells","ZR cells","HSP+ cells"),
                         colour = my.cell.type.cols))
temp_cellprop_df_all <- NULL
for(i in c(1:length(temp_Names))) {
  temp_filtered_summary <- subset(dataset@meta.data, orig.ident == temp_Names[i])
  temp_cellprop_df <- data.frame(unclass(table(temp_filtered_summary[,temp_colname])))
  temp_cellprop_df_all <- rbind(temp_cellprop_df_all,
                                (temp_cellprop_df <- data.frame(sample = rep(temp_Names[i], times = length(row.names(temp_cellprop_df))),
                                                                cell_type = row.names(temp_cellprop_df), value = temp_cellprop_df[,1],
                                                                proportions = temp_cellprop_df[,1]/colSums(temp_cellprop_df)))) }
temp_cellprop_df_all$cell_type <- factor(temp_cellprop_df_all$cell_type, levels=colour_pal$celltype)

prop_ggplot <- ggplot(temp_cellprop_df_all, aes(x=sample,fill=cell_type, y=proportions)) +
  geom_bar(stat="identity", color="black") +
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
        legend.text = element_text(size=7) + 
  xlab("Patient ID") +
  ylab("Cell Proportions")  +
  guides(fill=guide_legend(title="Cell Type")) +
  scale_fill_manual(values = as.vector(colour_pal$colour))


### Feature Plot for marker genes

DefaultAssay(dataset.combined)<- "RNA"
p_FeaturePlot <- FeaturePlot(dataset.combined, features = features, ncol = 4, cols=c("lightgrey", "gold","darkred"), min.cutoff = 'q5', max.cutoff = 'q95')  


### DotPlot Enrichment analyses

dotplot(Reactome_enricher, showCategory = 3, font.size = 10)


### Trajectory DDRTree

p_DDRT_Pseudotime <- plot_cell_trajectory(cds_fibro, color_by = "Pseudotime", show_branch_points = T, show_tree = TRUE, cell_size = 0.5) +  
        scale_color_viridis()

p_DDRT_Celltype <-  plot_cell_trajectory(cds_fibro, color_by = "cell_type", cell_size = 0.5) + 
  scale_color_manual(values = c("#28BBECFF", "#466BE3FF", "#30123BFF", "#FB8022FF", "#D23105FF", "#7A0403FF"))


### Trajectory Genes in Pseudotime

my_genes_1 <- c("NGF", "RGS5", "SEMA5A", "CD36")
my_genes_2 <- c("BNC2", "BICC1", "EFEMP1", "VCAN") 

cds_subset_1 <- cds_fibro[my_genes_1,]
p_genes_branch_1 <- plot_genes_branched_pseudotime(cds_subset_1, branch_point = 1, color_by = "cell_type", ncol = 4, panel_order = c("RGS5", "SEMA5A", "NGF", "CD36")  ) +
  scale_color_manual(values =  c("#28BBECFF", "#466BE3FF", "#30123BFF", "#FB8022FF", "#D23105FF", "#7A0403FF")) 

cds_subset_2 <- cds_fibro[my_genes_2,]
p_genes_branch_2 <- plot_genes_branched_pseudotime(cds_subset_2, branch_point = 1, color_by = "cell_type", ncol = 4, panel_order = c("BNC2", "BICC1", "EFEMP1", "VCAN")  ) +
  scale_color_manual(values =  c("#28BBECFF", "#466BE3FF", "#30123BFF", "#FB8022FF", "#D23105FF", "#7A0403FF")) 
  
              


