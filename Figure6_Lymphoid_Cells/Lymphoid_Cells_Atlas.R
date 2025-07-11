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


### Subset lymphoid cells

Idents(dataset) <- dataset$celltype_clusters
dataset_lympho <- subset(dataset, idents="Lymphoid cells")
Idents(dataset_lympho) <- dataset_lympho$celltype_ScibetGarnett
dataset_lympho <- subset(dataset_lympho, idents=c("B.cell","NK","T.cell")) 


### Normalization

DefaultAssay(dataset_lympho) <- "RNA"
dataset_lympho[['SCT']] <- NULL
dataset_lympho <- NormalizeData(dataset_lympho)
dataset_lympho <- FindVariableFeatures(dataset_lympho)
dataset_lympho <- ScaleData(dataset_lympho)
dataset_lympho <- SCTransform(dataset_lympho, vars.to.regress="percent.mt", verbose = TRUE)


### Dimensionality Reduction - Clustering

# PCA
dataset_lympho <- RunPCA(dataset_lympho, npcs = 30, verbose = TRUE)
ElbowPlot(dataset_lympho, ndims = 30)

# UMAP & Clustering
dataset_lympho <- RunUMAP(dataset_lympho, dims = 1:5, verbose = FALSE)
dataset_lympho <- FindNeighbors(dataset_lympho, dims = 1:5, verbose = FALSE)
dataset_lympho <- FindClusters(dataset_lympho, resolution = c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4))

# Test clustering stability
clustering_info <- dataset_lympho@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")
Idents(dataset_lympho) <- dataset_lympho$SCT_res.0.2

# UMAP visualization
DimPlot(dataset_lympho, reduction = "umap", label = TRUE)
DimPlot(dataset_lympho, reduction = "umap", group.by = "orig.ident", label = TRUE)


### Clusters identification

# DEG
markers_cluster <- FindAllMarkers(dataset_lympho, assay="RNA", slot="counts", test.use="wilcox", logfc.threshold=0.25, only.pos=T) # top100 in Suppl Table 12

# Rename clusters
new.cluster.ids <- c("Exhausted T cells", "Unassigned T cells", "Naive/Memory T cells", "NK-like T cells", "Plasma B cells")
names(new.cluster.ids) <- levels(dataset_lympho)
dataset_lympho <- RenameIdents(dataset_lympho, new.cluster.ids)
dataset_lympho$cell_type <- factor(dataset_lympho@active.ident)
saveRDS(dataset_lympho, "/path/to/Seurat_objects/Merge/dataset_lymphoid_cells.RDS")


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

saveRDS(Go_enricher, "/path/to/enrichment/results/GO_enricher_lymphoid_cells.RDS")
