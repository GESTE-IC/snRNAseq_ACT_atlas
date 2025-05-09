library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(cowplot)
library(future)
library(RColorBrewer)
library(tidyverse)
library(scibet)
library(viridis)
library(ggsci)
library(monocle)
library(garnett)
library(org.Hs.eg.db)
library(infercnv)
library(clustree)

dataset <- readRDS("/path/to/Seurat/objects/Merge/dataset_filtered.RDS")


### Normalization - Dimensionality Reduction - Clustering

# SCT Normalization
dataset <- SCTransform(dataset, conserve.memory = TRUE, verbose = TRUE)

# PCA
dataset <- RunPCA(dataset, npcs = 100, verbose = TRUE)
ElbowPlot(dataset, ndims = 100)

# UMAP & Clustering
dataset <- RunUMAP(dataset, dims = 1:50, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:50, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = c(0.8,1.0,1.2,1.4,1.6,1.8,2.0))
clustering_info <- dataset@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")

