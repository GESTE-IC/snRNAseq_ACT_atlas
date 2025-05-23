library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(future)
library(clustree)
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

image.dir = "/path/to/spaceranger_count/NAD4/outs/spatial/"
data.dir = "/path/to/spaceranger_count/NAD4/outs/"
NAD4 = Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "tissue_hires_image.png",
  filter.matrix = TRUE,
  to.upper = FALSE
)

image.dir = "/path/to/spaceranger_count/NAD5/outs/spatial/"
data.dir = "/path/to/spaceranger_count/NAD5/outs/"
NAD5 = Load10X_Spatial(
  data.dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "tissue_hires_image.png",
  filter.matrix = TRUE,
  to.upper = FALSE
)

NAD4$orig.ident <- "NAD4"
NAD5$orig.ident <- "NAD5"


### Integration

nadlist <- list(NAD4, NAD5)
names_nadlist <- c("NAD4", "NAD5")
nadlist <- setNames(nadlist, names_nadlist)

nadlist <- PrepSCTIntegration(object.list = nadlist, anchor.features = features)
nadlist <- lapply(X = nadlist, FUN = RunPCA, features = features)

anchors <- FindIntegrationAnchors(object.list = nadlist, normalization.method = "SCT", anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
nadintegrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)


### Dimensionality Reduction - Clustering

# PCA
nadintegrated <- RunPCA(nadintegrated, verbose = FALSE)

# UMAP & Clustering
nadintegrated <- RunUMAP(nadintegrated, reduction = "pca", dims = 1:30)
nadintegrated <- FindNeighbors(nadintegrated, reduction = "pca", dims = 1:30, verbose = FALSE)
for (i in res) {  nadintegrated <- FindClusters(object = nadintegrated, resolution = i) }

# Test clustering stability
clustree <- clustree(nadintegrated)
Idents(nadintegrated) <- nadintegrated$integrated_snn_res.0.2

# DEG - top 100 in Suppl Table 5
markersintegrated <- FindAllMarkers(nadintegrated, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = TRUE, assay = "Spatial", slot = "count")

nadintegrated[["old.ident.clusters"]] <- Idents(nadintegrated)
nadintegrated <- RenameIdents(nadintegrated, 
                              `0` = "ZF", 
                              `1` = "Capsule", 
                              `2` = "ZF HSP+",
                              `3` = "ZR", 
                              `4` = "ZG", 
                              `5` = "ZR / Medulla", 
                              `6` = "Medulla / ZR - Macrophages enriched")
nadintegrated$new.ident.clusters <- nadintegrated@active.ident
nadintegrated$new.ident.clusters <- factor(nadintegrated$new.ident.clusters, levels = c("Capsule", "ZG", "ZF", "ZF HSP+", "ZR", "ZR / Medulla", "Medulla / ZR - Macrophages enriched"))

saveRDS(nadintegrated, "/path/to/results/nadintegrated_annot.RDS")


### Add signatures scores

# deconvolution NAd4 - cell2location
cell2loc_ab <- read.table("/path/to/deconvolution/folder/results/cell2location_map/abundance_matrix.csv", sep = ",", header = TRUE)
cell2loc_ab <- cell2loc_ab %>%
  column_to_rownames("sample_name")
colnames(cell2loc_ab) <- substring(colnames(cell2loc_ab), 24)
cell2loc <- as.data.frame(t(apply(cell2loc_ab, 1, function(x) { x/sum(x) })))
colnames(cell2loc) <- paste0("prop.", colnames(cell2loc))
NAD4 <- AddMetaData(NAD4, metadata = cell2loc)
saveRDS(NAD4, "/path/to/cell2location_results/NAD4_deconv.RDS")

# deconvolution NAd5 - cell2location
cell2loc_ab <- read.table("/path/to/cell2location_results/cell2location_map/abundance_matrix.csv", sep = ",", header = TRUE)
cell2loc_ab <- cell2loc_ab %>%
  column_to_rownames("sample_name")
colnames(cell2loc_ab) <- substring(colnames(cell2loc_ab), 24)
cell2loc <- as.data.frame(t(apply(cell2loc_ab, 1, function(x) { x/sum(x) })))
colnames(cell2loc) <- paste0("prop.", colnames(cell2loc))
NAD5 <- AddMetaData(NAD5, metadata = cell2loc)
saveRDS(NAD5, "/path/to/cell2location_results/NAD5_deconv.RDS")


