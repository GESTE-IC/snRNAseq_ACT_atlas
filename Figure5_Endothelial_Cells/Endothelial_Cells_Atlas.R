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


### Subset endoth cells

Idents(dataset) <- dataset$celltype_clusters
dataset_endoth <- subset(dataset, idents="Endothelial cells") 
Idents(dataset_endoth) <- dataset_endoth$celltype_ScibetGarnett
dataset_endoth <- subset(dataset_endoth, idents="Endothelial.cell") 


### Normalization & Integration 

DefaultAssay(dataset_endoth) <- "RNA"
dataset_endoth[['SCT']] <- NULL
dataset_endoth <- NormalizeData(dataset_endoth)
dataset_endoth <- FindVariableFeatures(dataset_endoth)
dataset_endoth <- ScaleData(dataset_endoth)
dataset_endoth <- SCTransform(dataset_endoth, vars.to.regress="percent.mt", verbose = TRUE)


### Dimensionality Reduction - Clustering

# PCA
dataset_endoth <- RunPCA(dataset_endoth, npcs = 30, verbose = TRUE)
ElbowPlot(dataset_endoth, ndims = 30)

# UMAP & Clustering
dataset_endoth <- RunUMAP(dataset_endoth, dims = 1:10, verbose = FALSE)
dataset_endoth <- FindNeighbors(dataset_endoth, dims = 1:10, verbose = FALSE)
dataset_endoth <- FindClusters(dataset_endoth, resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8))

# Test clustering stability
clustering_info <- dataset_endoth@meta.data
clustree(clustering_info, prefix="SCT_snn_res.")
Idents(dataset_endoth) <- dataset_endoth$integrated_snn_res.0.1

# UMAP visualization
DimPlot(dataset_endoth, reduction = "umap", label = TRUE)
DimPlot(dataset_endoth, reduction = "umap", group.by = "orig.ident", label = TRUE)


### Clusters identification

# DEG
markers_cluster <- FindAllMarkers(dataset_endoth, assay="RNA", slot="counts", test.use="wilcox", logfc.threshold=0.25, only.pos=T) # top100 in Suppl Table 10

# Rename clusters
new.cluster.ids <- c("EC-venous","TEC1","TEC2","EC-HSP+","EC-arterial","EC-lymphatic")
names(new.cluster.ids) <- levels(dataset_endoth)
dataset_endoth <- RenameIdents(dataset_endoth, new.cluster.ids)
dataset_endoth$cell_type <- factor(dataset_endoth@active.ident, levels = c("EC-lymphatic","EC-arterial","EC-HSP+","EC-venous","TEC1","TEC2"))
saveRDS(dataset_endoth, "/path/to/Seurat_objects/Merge/dataset_endothelial_cells.RDS")


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

saveRDS(Go_enricher, "/path/to/enrichment/results/GO_enricher_endothelial_cells.RDS")
