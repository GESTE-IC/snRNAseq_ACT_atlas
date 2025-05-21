library(Matrix)
library(dplyr)
library(Seurat)
library(patchwork)
library(sctransform)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(clustree)

dataset <- readRDS("/path/to/Seurat_objects/Merge/dataset_filtered_normalized_annotated.RDS")


### Subset fibroblasts

Idents(dataset) <- dataset$celltype_clusters
dataset_fibro <- subset(dataset, idents="Fibroblasts") 
Idents(dataset_fibro) <- dataset_fibro$celltype_ScibetGarnett
dataset_fibro <- subset(dataset_fibro, idents="CAF") 


### Normalization

DefaultAssay(dataset_fibro) <- "RNA"
dataset_fibro[['SCT']] <- NULL
dataset_fibro <- NormalizeData(dataset_fibro)
dataset_fibro <- FindVariableFeatures(dataset_fibro)
dataset_fibro <- ScaleData(dataset_fibro)
dataset_fibro <- SCTransform(dataset_fibro, vars.to.regress="percent.mt", verbose = TRUE)


### Dimensionality Reduction - Clustering

# PCA
dataset_fibro <- RunPCA(dataset_fibro, npcs = 30, verbose = TRUE)
ElbowPlot(dataset_fibro, ndims = 30)

# UMAP & Clustering
dataset_fibro <- RunUMAP(dataset_fibro, dims = 1:10, verbose = FALSE)
dataset_fibro <- FindNeighbors(dataset_fibro, dims = 1:10, verbose = FALSE)
dataset_fibro <- FindClusters(dataset_fibro, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))

# Test clustering stability
clustering_info <- dataset_fibro@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")
Idents(dataset_fibro) <- dataset_fibro$integrated_snn_res.0.3

# UMAP visualization
DimPlot(dataset_fibro, reduction = "umap", label = TRUE)
DimPlot(dataset_fibro, reduction = "umap", group.by = "orig.ident", label = TRUE)


### Clusters identification

# DEG
markers_cluster <- FindAllMarkers(dataset_fibro, assay="RNA", slot="counts", test.use="wilcox", logfc.threshold=0.25, only.pos=T) 

# Rename clusters
new.cluster.ids <- c("Resident fibroblasts 1","CAF1","CAF3","CAF2",
                     "Resident fibroblasts 2","Resident fibroblasts 3",
                     "Pericytes", "Adipocytes", "Schwann cells")
names(new.cluster.ids) <- levels(dataset_fibro)
dataset_fibro <- RenameIdents(dataset_fibro, new.cluster.ids)
dataset_fibro$cell_type <- dataset_fibro@active.ident
saveRDS(dataset_fibro, "/path/to/Seurat_objects/Merge/dataset_fibroblasts.RDS")

# Aggregate clusters
new.cluster.ids <- c("Resident fibroblasts","CAF","CAF","CAF","Resident fibroblasts","Resident fibroblasts", "Pericytes", "Adipocytes", "Schwann cells")
names(new.cluster.ids) <- levels(dataset_fibro)
dataset_fibro <- RenameIdents(dataset_fibro, new.cluster.ids)

# DEG on aggregated clusters
markers_cluster <- FindAllMarkers(dataset_fibro, assay="RNA", slot="counts", test.use="wilcox", logfc.threshold=0.25, only.pos=T) # top100 in Suppl Table 8

# Overrepresentation analysis

markers_cluster <- markers_cluster[which(markers_cluster$p_val_adj<0.05),]
markers_cluster$cluster <- as.character(markers_cluster$cluster)

df_go <- as.data.frame(org.Hs.egGO)
go_universe <- unique(sort(df_go$gene_id))
GO_geneset <-  msigdbr(species = "Homo sapiens", category="C5", subcategory="GO:BP")
gs2gene_Go <- GO_geneset[, c("gs_name", "entrez_gene")]
gs2name_Go <- GO_geneset[, c("gs_name", "gene_symbol")]

all_clusters_markers <- NULL
for (cluster in unique(markers_cluster$cluster)) {
  name_cluster <- gsub(" ", "_", cluster)
  name_cluster <- gsub("-", "_", name_cluster)
  name_cluster <- gsub("\\+", "", name_cluster)
  list_markers <- DEG[which(DEG$cluster == cluster), "gene"][1:100]
  list_markers <- list_markers[which(!is.na(list_markers))]
  genes_EntrezID <- bitr(list_markers, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  cluster_markers <- data.frame(EntrezID = genes_EntrezID$ENTREZID, clusterID = cluster)
  all_clusters_markers <- rbind(all_clusters_markers, cluster_markers)   }

Go_enricher <- compareCluster(EntrezID ~ clusterID,  data = all_clusters_markers, fun = "enricher", TERM2GENE = gs2gene_Go, TERM2NAME = gs2name_Go)
Go_enricher@compareClusterResult$Description <- Go_enricher@compareClusterResult$ID

saveRDS(Go_enricher, "/path/to/enrichment/results/GO_enricher_fibroblasts.RDS")
