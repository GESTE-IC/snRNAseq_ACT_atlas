library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(viridis)
library(clustree)

dataset <- readRDS("/path/to/Seurat/objects/Merge/dataset_filtered_normalized_annotated.RDS")


### Violin Plots QC

my.cell.type.cols <- c("Myeloid cells" = "#21908CFF", "Lymphoid cells" = "#5DC863FF", "Steroid cells" = "#FDE725FF", 
                       "Chromaffin cells" = "violetred3",  "Endothelial cells" = "#3B528BFF", "Fibroblasts" = "#440154FF")

p_counts_samples <- VlnPlot(dataset, features = "nCount_RNA", group.by = "orig.ident", pt.size = 0, split.by = "histotype") +
  scale_fill_manual(values = c("darkmagenta", "hotpink3", "deepskyblue1", "mediumaquamarine", "gold"))

p_features_samples <- VlnPlot(dataset, features = "nFeature_RNA", group.by = "orig.ident", pt.size = 0, split.by = "histotype") +
  scale_fill_manual(values = c("darkmagenta","hotpink3","deepskyblue1","mediumaquamarine","gold"))

p_counts_celltype <- VlnPlot(dataset, features = "nCount_RNA", group.by = "celltype_clusters", cols = my.cell.type.cols, pt.size = 0)

p_features_celltype <- VlnPlot(dataset, features = "nFeature_RNA", group.by = "celltype_clusters", cols = my.cell.type.cols, pt.size = 0)


### Clustree

clustering_info <- dataset@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")


### UMAP

p_umap_clusters <- DimPlot(dataset, reduction = "umap", group.by = "SCT_snn_res.1.2", label = TRUE, repel = FALSE, raster = FALSE) + 
  scale_color_viridis_d(option = "turbo")

p_umap_celltype <- DimPlot(dataset, reduction = "umap", group.by = "celltype_clusters", label = FALSE, raster = FALSE, cols = my.cell.type.cols)

p_umap_sample <- DimPlot(dataset, reduction = "umap", group.by = "orig.ident", raster = FALSE, label = FALSE) + 
  scale_color_viridis_d(option = "plasma") +
  ggtitle("")

p_umap_histotype <- DimPlot(dataset, reduction = "umap", group.by = "histotype", raster = FALSE, label = FALSE,
              cols = c("darkmagenta", "hotpink3", "deepskyblue1", "mediumaquamarine", "gold")) +
  ggtitle("")

p_umap_CNA <- DimPlot(dataset, reduction = "umap", group.by = "CNA_status_steroid", raster=F, label=F,
             cols=c("gold", "mediumpurple4"), na.value = "grey70") +
  ggtitle("")

p_umap_scibet_garnett <- DimPlot(dataset, reduction = "umap", group.by = "celltype_ScibetGarnett", raster=F, label=F) +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  ggtitle("")


### Stacked Violin Plot for marker genes

features <- c("SCARB1", "CYP11A1", "STAR","NR5A1",
               "TH","DBH","CHGA",
               "FLT1","VWF","PTPRB", 
               "LAMA2","FN1","SERPINE1",
               "CD163","MSR1","F13A1",
               "IKZF1","CD247","THEMIS")

dataset$celltype_reorder <- factor(dataset$celltype_clusters, levels = c("Lymphoid cells", "Myeloid cells",  "Fibroblasts", "Endothelial cells", "Chromaffin cells", "Steroid cells"))

p_Vln_stacked <- VlnPlot(dataset, features = features,  cols = my.cell.type.cols, fill.by = 'ident', assay = "RNA", pt.size =, 0, adjust=1.5, stack = TRUE) +
  ggtitle("") +
  scale_y_discrete(limits = levels(dataset$celltype_reorder)) +
  NoLegend() +
  theme(plot.title = element_text(size = 15, hjust = 0.5), axis.title = element_blank())


### Proportion of cells in each sample

temp_Names <- unique(dataset$orig.ident)
dataset@meta.data$cell_type <- dataset@active.ident
temp_colname <- "cell_type"
colour_pal <- data.frame(celltype = c("Fibroblasts","Endothelial cells","Steroid cells", "Chromaffin cells","Myeloid cells","Lymphoid cells"),
                                colour = c("#FDE725FF","violetred3","#3B528BFF","#440154FF","#21908CFF","#5DC863FF"))
temp_cellprop_df_all <- NULL
for(i in c(1:length(temp_Names))) {
  temp_filtered_summary <- subset(dataset@meta.data, orig.ident == temp_Names[i])
  temp_cellprop_df <- data.frame(unclass(table(temp_filtered_summary[,temp_colname])))
  temp_cellprop_df_all <- rbind(temp_cellprop_df_all,
                                (temp_cellprop_df <- data.frame(sample = rep(temp_Names[i], times = length(row.names(temp_cellprop_df))),
                                                                cell_type = row.names(temp_cellprop_df), value = temp_cellprop_df[,1],
                                                                proportions = temp_cellprop_df[,1]/colSums(temp_cellprop_df),
                                                                subtype = unique(temp_filtered_summary$histotype)))) }
temp_cellprop_df_all$cell_type <- factor(temp_cellprop_df_all$cell_type, levels=colour_pal$celltype)
temp_cellprop_df_all$subtype <- factor(temp_cellprop_df_all$subtype, levels=c("ACC C1A","ACC C1B","ACA", "PBMAH","Normal adrenal"))

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
  scale_fill_manual(values = as.vector(colour_pal$colour)) +
  facet_grid(. ~ subtype, scales="free", space = "free")

