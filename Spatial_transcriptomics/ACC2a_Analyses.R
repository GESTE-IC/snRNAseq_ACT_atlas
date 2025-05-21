library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(tibble)
library(ggpubr)
library(clustree)


### Create Seurat object

image.dir = "/path/to/spaceranger_count/ACC2a/outs/spatial/"
data.dir = "/path/to/spaceranger_count/ACC2a/outs/"
ACC2a = Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "tissue_hires_image.png",
  filter.matrix = TRUE,
  to.upper = FALSE
)

figure_dir <- "/path/to/figures/"
sample_name <- "ACC2a_"

ACC2a$orig.ident = "ACC2a"
ACC2a$techno.ident <- "Spatial_Visium"


### Normalization 

ACC2a[["percent.mt"]] <- PercentageFeatureSet(ACC2a, pattern = "^MT-")
ACC2a <- SCTransform(ACC2a, assay = "Spatial", verbose = TRUE, return.only.var.genes = FALSE, vars.to.regress = "percent.mt")


### Dimensionality Reduction - Clustering

# PCA
ACC2a <- RunPCA(ACC2a, assay = "SCT", npcs = 50, verbose = TRUE, features = VariableFeatures(ACC2a))

# UMAP & Clustering
ACC2a <- RunUMAP(ACC2a, reduction = "pca", dims = 1:30, verbose = TRUE)
ACC2a <- FindNeighbors(ACC2a, reduction = "pca", dims = 1:30, verbose = TRUE)
res <- seq(0.1, 1, 0.1)
for (i in res) {
  ACC2a <- FindClusters(object = ACC2a, resolution = i)
}

# Test clustering stability
clustree <- clustree(ACC2a)
Idents(ACC2a) <- ACC2a$SCT_snn_res.0.3

new.cluster.ids <- c("cluster1", "cluster2", "cluster3")
names(new.cluster.ids) <- levels(ACC2a)
ACC2a <- RenameIdents(ACC2a, new.cluster.ids)

# DEG - top 100 in Suppl Table 17
ACC2a_markers <- FindAllMarkers(ACC2a, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = TRUE, assay = "Spatial", slot = "data")


### Add signatures scores

# deconvolution TMA - cell2location
cell2loc_ab <- read.table("/path/to/cell2location_results/cell2location_map/abundance_matrix.csv", sep = ",", header = TRUE)
cell2loc_ab <- cell2loc_ab %>%
  column_to_rownames("sample_name")
colnames(cell2loc_ab) <- substring(colnames(cell2loc_ab), 24)
cell2loc <- as.data.frame(t(apply(cell2loc_ab, 1, function(x) { x/sum(x) } )))
ACC2a <- AddMetaData(ACC2a, metadata = cell2loc)

# gene modules scores
gene_modules <- readRDS("/path/to/GeneModules/results/Gene_modules.RDS"))
gene_modules <- gene_modules[order(names(gene_modules))]
ACC2a <- AddModuleScore(ACC2a, gene_modules, name = paste0(names(gene_modules), "_"))

saveRDS(ACC2a, "/path/to/results/ACC2a_annot.RDS")
