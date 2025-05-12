library(Matrix)
library(dplyr)
library(Seurat)
#library(SeuratData)
library(patchwork)
library(ggplot2)
library(monocle)
library(SeuratWrappers)

datase_fibro <- readRDS("/path/to/Seurat_objects/Merge/dataset_fibroblasts.RDS")
dataset_traj <- subset(dataset_fibro, idents = c("Pericytes","Adipocytes","Schwann cells"), invert=T)


### Ordering algorithm using variable genes

# select genes
ordering_genes <- VariableFeatures(dataset_traj)[1:1000]

# Create cds
data <- as(as.matrix(dataset_traj@assays$SCT@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = dataset_traj@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = uninormal())
cds <- setOrderingFilter(cds, ordering_genes)

# reduce dimensionality
cds <- reduceDimension(cds, norm_method="none", reduction_method="DDRTree")

# trajectory
cds <- orderCells(cds)

saveRDS(cds, "/path/to/Monocle/results/Cds_Fibroblasts.RDS")


### DEG with pseudotime

diff_test_res_all <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 1)
diff_test_res_all <- diff_test_res_all[order(diff_test_res_all$qval),] # Suppl Table 9

saveRDS(diff_test_res_all, "/path/to/Monocle/results/DEG_byPseudotime_Fibroblasts.RDS")
