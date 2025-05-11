library(Matrix)
library(dplyr)
library(Seurat)

dataset <- readRDS("/path/to/Seurat/objects/Merge/dataset.RDS")

# Create matrices for Scrublet

path_input_scrublet <- "/path/to/Scrublet/input/"

for (sample in unique(dataset$orig.ident)) {
  IDs_current_sample <- row.names(dataset@meta.data[which(dataset$orig.ident==sample), ])
  writeMM(dataset@assays$RNA@counts[, IDs_current_sample], paste0(path_input_scrublet, "Rawdata_", sample, "_", date, ".mtx"))
}

# Remove doublets using Scrublet output

path_output_scrublet <- "/path/to/Scrublet/outputs/"
list_outputs <- list.files(path_output_scrublet, pattern=NULL, all.files=FALSE, full.names=FALSE)

list_singletons <- NULL
for (output in list_outputs) {
  sample_name <- gsub("Rawdata_", "", output)
  sample_name <- gsub(paste0("_table.csv"), "", sample_name)
  print(sample_name)
  ID_cells_sample <- row.names(dataset@meta.data[which(dataset@meta.data$orig.ident==sample_name), ])
  sample_score <- read.csv(file = paste0(path_output_scrublet,"/",output), header=TRUE)
  sample_singletons <- ID_cells_sample[which(sample_score$predicted_doublet=="False")]
  list_singletons <- c(list_singletons, sample_singletons)
}

dataset <- dataset[,list_singletons]

# Filtering on nFeature & percent.mt

dataset[["percent.mt"]] <- PercentageFeatureSet(dataset, pattern = "^MT-")
dataset <- subset(dataset, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 5)

saveRDS(dataset, "/path/to/Seurat/objects/Merge/dataset_filtered.RDS")
