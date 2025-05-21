library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(RColorBrewer)
library(pheatmap)

### Add TME annotations in whole dataset

dataset <- readRDS("/path/to/Seurat_objects/Merge/dataset_filtered_normalized_annotated_GM.RDS")

dataset_fibro <- readRDS("/path/to/Seurat_objects/Merge/dataset_fibroblasts.RDS")
dataset_endoth <- readRDS("/path/to/Seurat_objects/Merge/dataset_endothelial_cells.RDS")
dataset_lympho <- readRDS("/path/to/Seurat_objects/Merge/dataset_lymphoid_cells.RDS")
dataset_myelo <- readRDS("/path/to/Seurat_objects/Merge/dataset_myeloid_cells.RDS")

dataset$celltype_tme <- as.character(dataset@active.ident)
dataset$celltype_tme[which(dataset$celltype_tme != "Steroid cells")] <- NA
dataset$celltype_tme[which(dataset$celltype_ScibetGarnett != "Steroid cells")] <- NA
dataset$celltype_tme[names(dataset_endoth$cell_type)] <- as.character(dataset_endoth$cell_type)
dataset$celltype_tme[names(dataset_fibro$cell_type)] <- as.character(dataset_fibro$cell_type)
dataset$celltype_tme[names(dataset_myelo$cell_type)] <- as.character(dataset_myelo$cell_type)
dataset$celltype_tme[names(dataset_lympho$cell_type)] <- as.character(dataset_lympho$cell_type)


### Select annotations for deconvolution 
# remove cell type with too few cells ()
# keep only begining and end of trajectories

dataset$celltype_deconv <- dataset$celltype_tme
dataset$celltype_deconv[which(dataset$celltype_deconv %in% c("EC-lymphatic", "EC-HSP+", "TEC1",
  "CAF2", "Resident fibroblasts 2","Resident fibroblasts 3", "Pericytes", "Adipocytes", "Schwann cells",
  "Resident macrophages 2", "TAM1", "Cycling macrophages", "Unassigned T cells", "Plasma B cells") ] <- NA

saveRDS(dataset, "/path/to/Seurat_objects/Merge/dataset_filtered_normalized_annotated_GM_TME.RDS")
                              
                              
### Split dataset for deconvolution and prepare matrix for cibersort

# counts table
counts<- dataset@assays$RNA@counts

# selection of celltypes for deconvolution
dataset@meta.data$celltype_deconv <- gsub(" " , "_", dataset@meta.data$celltype_deconv)
dataset@meta.data$celltype_deconv <- gsub("/" , "_", dataset@meta.data$celltype_deconv)
dataset@meta.data$celltype_deconv <- gsub("-" , "_", dataset@meta.data$celltype_deconv)
celltypes <- unique(dataset@meta.data$celltype_deconv)

# random split                         
cells_train <- NULL
cells_test <- NULL
for (celltype in celltypes) {
  cells <- rownames(dataset@meta.data)[which(dataset@meta.data$celltype_deconv == celltype)]
  cells1 <- sample(cells, min(1000,floor(length(cells)/2)))
  cells2 <- setdiff(cells, cells1)
  cells_train <- c(cells_train, cells1)
  cells_test <- c(cells_test, cells2)
}

# create training dataset
dataset_train <- subset(dataset, cells=cells_train)
counts_train <- as.matrix(dataset_train@assays$RNA@data)
counts_train <- cbind(rownames(counts_train),counts_train)
colnames(counts_train) <- c("Gene", dataset_train@meta.data[,"celltype_deconv"])
write.table(counts_train, file = "path/to/cibersort/ref_matrices_pseudobulk/refmatrix_celltype_deconv.txt", sep="\t", dec =".", col.names = T, row.names = F, quote = F)

# create pseudo-bulk validation dataset
dataset_test <- subset(dataset, cells=cells_test)
real_prop <- prop.table(table(as.character(dataset_test$orig.ident), dataset_test$celltype_deconv),1)
write.table(real_prop, "path/to/cibersort/results/Real_prop_celltype_deconv.txt")
counts_test <- AverageExpression(dataset_test, assays = "RNA", group.by = "orig.ident", slot="data")
counts_test <- counts_test$RNA
colnames(counts_test) <- unname(colnames(counts_test))
samples <- colnames(counts_test)
counts_test <- as.matrix(counts_test)
counts_test<- cbind(rownames(counts_test),counts_test)
colnames(counts_test)<- c("Gene", samples)
write.table(counts_test, file = "path/to/cibersort/ref_matrices_pseudobulk/pseudobulk_celltype_deconv.txt", sep="\t", dec =".", col.names = T, row.names = F, quote = F)


### Pseudobulk results analyses

col <- colorRampPalette(c("blue", "white", "darkred"))(20)

# load results estimated proportions
res_cibersort <- read.csv2("path/to/cibersort/results/Pseudobulk_celltype_deconv/CIBERSORTx_Adjusted.txt", header=T, sep = "\t", dec=".")
rownames(res_cibersort)<- res_cibersort$Mixture
res_cibersort<- res_cibersort[,2:(ncol(res_cibersort)-3)]

# load real proportions
res_ref <- read.table("path/to/cibersort/results/Real_prop_celltype_deconv.txt") header=T)
res_ref_order<- res_ref[rownames(res_cibersort), colnames(res_cibersort)]

# rename cell pop
colnames(res_cibersort) <- paste0(colnames(res_cibersort),"_estimated")
colnames(res_ref_order) <- paste0(colnames(res_ref_order),"_observed")

# correlation plot
M_cor <- cor(res_cibersort, res_ref_order)
r_cor_global <- cor.test(as.matrix(res_cibersort), as.matrix(res_ref_order))
corrplot(M_cor, method="circle", order="hclust", type="upper", 
           tl.col="black", tl.srt=45, col=col,
           # title = paste("r =", round(r_cor_global$estimate,2)),
           mar = c(1, 1, 2, 1))
                            
