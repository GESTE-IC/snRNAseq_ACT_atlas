library(Matrix)
library(dplyr)
library(Seurat)
#library(SeuratData)
library(patchwork)
library(sctransform)
library(ggplot2)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(msigdbr)
library(clustree)

dataset <- readRDS("/path/to/Seurat_objects/Merge/dataset_filtered_normalized_annotated.RDS")


### Subset normal steroid cells

Idents(dataset) <- dataset$celltype_clusters
dataset_normaladrenal <- subset(dataset, idents="Steroid cells")
Idents(dataset_normaladrenal) <- dataset_normaladrenal$orig.ident
dataset_normaladrenal <- subset(dataset_normaladrenal, idents=c("NAd1", "NAd2", "NAd3", "NAd4"))
Idents(dataset_normaladrenal) <- dataset_normaladrenal$celltype_ScibetGarnett
dataset_normaladrenal <- subset(dataset_normaladrenal, idents="Steroid cells")


### Normalization & Integration 

# Normalize each sample independently with regress on percent.mt (otherwise associated with clustering)
dataset.list <- SplitObject(dataset_normaladrenal, split.by = "orig.ident")
dataset.list <- lapply(X = dataset.list, FUN = function(x) { x <- SCTransform(x, vars.to.regress = c("percent.mt"))})

# Integration using cca
features <- SelectIntegrationFeatures(object.list = dataset.list, nfeatures = 3000)
dataset.list <- PrepSCTIntegration(object.list = dataset.list, anchor.features = features)
dataset.anchors <- FindIntegrationAnchors(object.list = dataset.list, anchor.features = features, normalization.method = "SCT", reduction = "cca")
dataset.combined <- IntegrateData(anchorset = dataset.anchors, normalization.method = "SCT")
DefaultAssay(dataset.combined) <- "integrated"


### Dimensionality Reduction - Clustering

dataset.combined <- ScaleData(dataset.combined, verbose = FALSE, features= rownames(dataset.combined))

# PCA
dataset.combined <- RunPCA(dataset.combined, npcs = 50, verbose = FALSE)
ElbowPlot(dataset.combined, ndims=50)

# UMAP & Clustering
dataset.combined <- RunUMAP(dataset.combined, reduction = "pca", dims = 1:10)
dataset.combined <- FindNeighbors(dataset.combined, reduction = "pca", dims = 1:10)
dataset.combined <- FindClusters(dataset.combined, resolution = c(0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8))

# Test clustering stability
clustering_info <- dataset.combined@meta.data
clustree(clustering_info, prefix="integrated_snn_res.")
Idents(dataset.combined) <- dataset$integrated_snn_res.0.3

# UMAP visualization
DimPlot(dataset.combined, reduction = "umap", label = TRUE)
DimPlot(dataset.combined, reduction = "umap", group.by = "orig.ident", label = TRUE)


### Clusters identification

# DEG
markers_cluster <- FindAllMarkers(dataset.combined, assay="RNA", slot="counts", test.use="wilcox", logfc.threshold=0.25, only.pos=T) # top100 in Suppl Table 3

# Rename clusters
new.cluster.ids <- c("ZF cells","Intermediate-state cells","ZR cells", "ZG cells", "HSP+ cells", "ZG cells - CYP11B2+")
names(new.cluster.ids) <- levels(dataset.combined)
dataset.combined <- RenameIdents(dataset.combined, new.cluster.ids)
dataset.combined$celltype_integrated <- dataset.combined@active.ident
saveRDS(dataset.combined, "/path/to/Seurat_objects/Merge/dataset_normal_adrenal_integrated.RDS")

# Overrepresentation analysis

markers_cluster <- markers_cluster[which(markers_cluster$p_val_adj<0.05),]
markers_cluster$cluster <- as.character(markers_cluster$cluster)

Reactome_geneset <- msigdbr(species = "Homo sapiens", category="C2", subcategory="CP:REACTOME")
gs2gene_Reactome <- Reactome_geneset[, c("gs_name", "entrez_gene")]
gs2name_Reactome <- Reactome_geneset[, c("gs_name", "gene_symbol")]

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

Reactome_enricher <- compareCluster(EntrezID ~ clusterID,  data = all_clusters_markers, fun = "enricher", TERM2GENE = gs2gene_Reactome, TERM2NAME = gs2name_Reactome)
Reactome_enricher@compareClusterResult$Description <- Reactome_enricher@compareClusterResult$ID

saveRDS(Reactome_enricher, "/path/to/enrichment/results/Reactome_enricher_normal_adrenal.RDS")


### Normalization without Integration

dataset_normaladrenal$celltype_integrated <- dataset.combined@active.ident
DefaultAssay(dataset_normaladrenal) <- "RNA"
dataset_normaladrenal[['SCT']] <- NULL

# Normalization SCT
dataset_normaladrenal <- SCTransform(dataset_normaladrenal, verbose = FALSE)

# PCA
dataset_normaladrenal <- RunPCA(dataset_normaladrenal, npcs = 30, verbose = TRUE)
ElbowPlot(dataset_normaladrenal, ndims = 30)

# UMAP & Clustering
dataset_normaladrenal <- RunUMAP(dataset_normaladrenal, dims = 1:10, verbose = FALSE)
dataset_normaladrenal <- FindNeighbors(dataset_normaladrenal, dims = 1:10, verbose = FALSE)
dataset_normaladrenal <- FindClusters(dataset_normaladrenal, resolution = 0.8, verbose = FALSE)

# UMAP visualization
DimPlot(dataset_normaladrenal, reduction = "umap", label = TRUE) + NoLegend()
DimPlot(dataset_normaladrenal, reduction = "umap", group.by="orig.ident", label = TRUE)
DimPlot(dataset_normaladrenal, reduction = "umap", group.by="celltype_integrated", label = F)

saveRDS(dataset_normaladrenal, "/path/to/Seurat_objects/Merge/dataset_normal_adrenal_not_integrated.RDS"))
