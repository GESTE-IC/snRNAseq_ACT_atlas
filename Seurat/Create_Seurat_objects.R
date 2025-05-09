# Seurat v4.3.0

library(Matrix)
library(dplyr)
library(Seurat)

# create individual Seurat objects

input_path <- "/path/to/cellranger_output/"
files <- list.files(path=input_path)
output_path <- "/path/to/Seurat/objects/"

for (f in files) {
mat <- Read10X(data.dir=paste(input_path, f, "/outs/filtered_feature_bc_matrix/",sep=""))
obj <- CreateSeuratObject(counts=mat, project=f, min.cells=3, min.features=200)
saveRDS(obj,paste0(output_path,"Seurat_object_", f, ".RDS")) }

# merge Seurat objects into Seurat dataset

seurat_obj <- list.files(output_path)
names_seurat <- gsub("Seurat_object_", "", seurat_obj)
names_seurat <- gsub(".RDS", "", names_seurat)

obj1 <- readRDS(paste0(output_path, seurat_obj[1]))
list_to_merge <- vector("list", length(seurat_obj)-1)

for (i in 2:length(seurat_obj)){
  obj <- readRDS(paste0(output_path, seurat_obj[i]))
  list_to_merge[[i-1]] <- obj }

dataset <- merge(obj1, y=list_to_merge, add.cell.ids=c(names_seurat))
