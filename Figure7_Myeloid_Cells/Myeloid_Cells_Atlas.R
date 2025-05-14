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


### Subset myeloid cells

Idents(dataset) <- dataset$celltype_clusters
dataset_myelo <- subset(dataset, idents="Myeloid cells")
Idents(dataset_myelo) <- dataset_myelo$celltype_ScibetGarnett
dataset_myelo <- subset(dataset_myelo, idents="Macrophage")


### Normalization & Integration 

DefaultAssay(dataset_myelo) <- "RNA"
dataset_myelo[['SCT']] <- NULL
dataset_myelo <- NormalizeData(dataset_myelo)
dataset_myelo <- FindVariableFeatures(dataset_myelo)
dataset_myelo <- ScaleData(dataset_myelo)
dataset_myelo <- SCTransform(dataset_myelo, verbose = TRUE)


### Dimensionality Reduction - Clustering

# PCA
dataset_myelo <- RunPCA(dataset_myelo, npcs = 30, verbose = TRUE)
ElbowPlot(dataset_myelo, ndims = 30)

# UMAP & Clustering
dataset_myelo <- RunUMAP(dataset_myelo, dims = 1:10, verbose = FALSE)
dataset_myelo <- FindNeighbors(dataset_myelo, dims = 1:10, verbose = FALSE)
dataset_myelo <- FindClusters(dataset_myelo, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))

# Test clustering stability
clustering_info <- dataset_myelo@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")
Idents(dataset_myelo) <- dataset_myelo$integrated_snn_res.0.1

# UMAP visualization
DimPlot(dataset_myelo, reduction = "umap", label = TRUE)
DimPlot(dataset_myelo, reduction = "umap", group.by = "orig.ident", label = TRUE)


### Clusters identification

# DEG
markers_cluster <- FindAllMarkers(dataset_myelo, assay="RNA", slot="counts", test.use="wilcox", logfc.threshold=0.25, only.pos=T) # top100 in Suppl Table 14

# Rename clusters
new.cluster.ids <- c("Resident macrophages 1", "Inflammatory macrophages", "TAM1", "TAM2", "Resident macrophages 2","Perivascular TAM", "Cycling macrophages")
names(new.cluster.ids) <- levels(dataset_myelo)
dataset_myelo <- RenameIdents(dataset_myelo, new.cluster.ids)
dataset_myelo$cell_type <- factor(dataset_myelo@active.ident,  
                                  levels = c("Resident macrophages 1", "Resident macrophages 2", "Inflammatory macrophages", "TAM1", "TAM2","Perivascular TAM", "Cycling macrophages")
saveRDS(dataset_myelo, "/path/to/Seurat_objects/Merge/dataset_myeloid_cells.RDS")


# Overrepresentation analysis

dataset_myelo$ident_enrich[which(dataset_myelo$ident_enrich %in% c("TAM1", "TAM2"))] <- "TAM"
dataset_myelo$ident_enrich[which(dataset_myelo$ident_enrich %in% c("Resident macrophages 1", "Resident macrophages 2"))] <- "Resident macrophages"
Idents(dataset_myelo) <- dataset_myelo$ident_enrich
markers_cluster <- FindAllMarkers(dataset_myelo, only.pos = T, logfc.threshold = 0.25, min.pct = 0.1, test.use = "wilcox", assay="RNA", slot="counts")
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

saveRDS(Go_enricher, "/path/to/enrichment/results/GO_enricher_myeloid_cells.RDS")
