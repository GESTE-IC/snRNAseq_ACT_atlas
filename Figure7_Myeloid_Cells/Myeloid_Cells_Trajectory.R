library(Matrix)
library(dplyr)
library(Seurat)
library(SeuratData)
library(patchwork)
library(ggplot2)
library(monocle)
library(SeuratWrappers)

dataset_myelo <- readRDS("/path/to/Seurat_objects/Merge/dataset_myeloid_cells.RDS")
dataset_traj <- dataset_myelo


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

# Reorder with EC-venous as starting point
GM_state <- function(cds_data){ T0_counts <- table(pData(cds_data)$State, pData(cds_data)$cell_type)[,"Resident macrophages 1"]
return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))])) }
cds <- orderCells(cds, root_state = GM_state(cds))

saveRDS(cds, "/path/to/Monocle/results/Cds_Myeloid_Cells.RDS")


### DEG with pseudotime

BEAM_res_branch_1 <- BEAM(cds, branch_point = 2)

diff_test_res_all <- differentialGeneTest(cds, fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 1)
diff_test_res_all <- diff_test_res_all[order(diff_test_res_all$qval),] # Suppl Table 15

saveRDS(diff_test_res_all, "/path/to/Monocle/results/DEG_byPseudotime_Myeloid_Cells.RDS")
