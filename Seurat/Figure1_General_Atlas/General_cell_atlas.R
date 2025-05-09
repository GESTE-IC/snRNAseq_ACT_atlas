library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(tidyverse)
library(scibet)
library(ggsci)
library(monocle)
library(garnett)
library(org.Hs.eg.db)
library(clustree)

dataset <- readRDS("/path/to/Seurat/objects/Merge/dataset_filtered.RDS")


### Normalization - Dimensionality Reduction - Clustering

# Scale RNA assay for plots
dataset <- NormalizeData(dataset, normalization.method = "LogNormalize", scale.factor = 10000)

# SCT Normalization
dataset <- SCTransform(dataset, conserve.memory = TRUE, verbose = TRUE)

# PCA
dataset <- RunPCA(dataset, npcs = 100, verbose = TRUE)
ElbowPlot(dataset, ndims = 100)

# UMAP & Clustering
dataset <- RunUMAP(dataset, dims = 1:50, verbose = FALSE)
dataset <- FindNeighbors(dataset, dims = 1:50, verbose = FALSE)
dataset <- FindClusters(dataset, resolution = c(0.8,1.0,1.2,1.4,1.6,1.8,2.0))

# Test clustering stability
clustering_info <- dataset@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")
Idents(dataset) <- dataset$SCT_snn_res.1.2

# Add annotations
dataset$orig.ident <- factor(dataset$orig.ident, 
                  levels = c("ACC1a","ACC1b","ACC2a","ACC2b","ACC3","ACC4","ACC5","ACC6","ACC7","ACC8a","ACC8b",
                  "ACC9","ACC10","ACC11","ACC12","ACC13","ACC14","ACC15","ACC16","ACC17",
                  "ACA1","ACA2","ACA3","ACA4","ACA5","ACA6","ACA7","ACA8",
                  "PBMAH1","PBMAH2","PBMAH3","PBMAH4","PBMAH5","PBMAH6",
                  "NAd1","NAd2","NAd3","NAd4"))
dataset$histotype <- dataset$orig.ident 
levels(dataset$histotype) <- c(rep("ACC C1A",11), rep("ACC C1B",9), rep("ACA",8), rep("PBMAH",6), rep("Normal adrenal",4))
dataset$histotype <- factor(as.character(dataset$histotype), 
                    levels=c("ACC C1A", "ACC C1B", "ACA", "PBMAH", "Normal adrenal"))

# UMAP visualization
DimPlot(dataset, reduction = "umap", label = TRUE)
DimPlot(dataset, reduction = "umap", group.by = "orig.ident", label = T)
DimPlot(dataset, reduction = "umap", group.by = "histotype", label = FALSE)


### Clusters identification

# DEG
markers_cluster <- FindAllMarkers(dataset, only.pos=T, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox", assay="RNA", slot="counts")

# Garnett prediction for steroid and chromaffin cells
classifier <- readRDS("adrenal_classifier.RDS") # reference to be downloaded on https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/
data <- as(dataset@assays$RNA@data, 'dgTMatrix')
pd <- new('AnnotatedDataFrame', data = dataset@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.5,
                              expressionFamily = negbinomial.size())
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- classify_cells(monocle_cds, classifier,
                              db = org.Hs.eg.db,
                              cluster_extend = TRUE,
                              cds_gene_id_type = "SYMBOL")
cell_ID_garnett <- monocle_cds@phenoData@data$cluster_ext_type
names(cell_ID_garnett) <- rownames(monocle_cds@phenoData@data)
dataset$celltype_Scibet <- cell_ID_scibet[rownames(dataset@meta.data)]

# Scibet prediction for microenvironement cells
model <- readr::read_csv("GSE115978_normal_scibet_core.csv")
model<- pro.core(model)
genessel <- intersect(rownames(model),dataset@assays$RNA@data@Dimnames[[1]])
query <- t(as.matrix(dataset@assays$RNA@data[genessel,]))
prd <- LoadModel(model)
cell_ID_scibet <- prd(query)
names(cell_ID_scibet) <- rownames(query)
dataset$celltype_Garnett <- cell_ID_garnett[rownames(dataset@meta.data)]

# Garnett and Scibet
new_cell_ID <- cell_ID_scibet[rownames(dataset@meta.data)]
new_cell_ID[which(new_cell_ID=="T.cell"|new_cell_ID=="T.CD8"|new_cell_ID=="T.CD4")] <- "T.cell"
medulla_cell <- names(cell_ID_garnett)[which(cell_ID_garnett=="Chromaffin cells")]
new_cell_ID[medulla_cell] <- "Chromaffin cells"
steroid_cell <- names(cell_ID_garnett)[which(cell_ID_garnett=="Adrenocortical cells")]
new_cell_ID[steroid_cell] <- "Steroid cells"
dataset$celltype_ScibetGarnett <- new_cell_ID

# rename clusters by cell type
new.cluster.ids <- c("Myeloid cells","Endothelial cells","Steroid cells","Steroid cells","Steroid cells","Fibroblasts",
                     "Steroid cells","Steroid cells","Steroid cells","Steroid cells","Steroid cells",
                     "Steroid cells","Steroid cells","Steroid cells","Steroid cells","Steroid cells",
                     "Steroid cells","Steroid cells","Steroid cells","Steroid cells","Steroid cells",
                     "Steroid cells","Steroid cells","Steroid cells","Steroid cells","Myeloid cells",
                     "Myeloid cells","Steroid cells","Steroid cells","Steroid cells","Steroid cells",
                     "Steroid cells","Steroid cells","Steroid cells","Myeloid cells","Lymphoid cells",
                     "Steroid cells","Steroid cells","Steroid cells","Steroid cells","Steroid cells",
                     "Steroid cells","Fibroblasts","Steroid cells","Steroid cells","Steroid cells",
                     "Steroid cells","Endothelial cells","Steroid cells","Steroid cells","Steroid cells",
                     "Steroid cells","Fibroblasts","Myeloid cells","Fibroblasts","Steroid cells",
                     "Endothelial cells","Chromaffin cells","Steroid cells","Myeloid cells","Myeloid cells")
names(new.cluster.ids) <- levels(dataset)
dataset <- RenameIdents(dataset, new.cluster.ids)
dataset$celltype_clusters <- dataset@active.ident

# DEG for each cell type
markers_cluster <- FindAllMarkers(dataset, only.pos=T,  logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox", assay="RNA", slot="counts") # top 100 in Suppl Table 2

saveRDS(dataset, "/path/to/Seurat/objects/Merge/dataset_filtered_normalized.RDS")
