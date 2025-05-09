library(Matrix)
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(infercnv)

dataset <- readRDS("/path/to/Seurat/objects/Merge/dataset_filtered_normalized.RDS")


### Create gene ref file

gene_ordering_file <- read.table("/path/to/gencode_v21_gen_pos.complete.txt"))
genename <- strsplit(gene_ordering_file$V1,"[|]")
genename2 <- sapply(genename,function(x) x[1])
gene_ordering_file$V1 <- genename2

# remove genes not in dataset                  
selection <- intersect(gene_ordering_file$V1,rownames(dataset@assays$RNA@counts))
gene_ordering_file <- gene_ordering_file[which(gene_ordering_file$V1 %in% selection),]

# remove duplicates                    
gene_ordering_file <- gene_ordering_file[which(!duplicated(gene_ordering_file[,1])),] 

# order by chr pos                    
gene_ordering_file <- gene_ordering_file[order(gene_ordering_file$V3),] 
gene_ordering_file$V2 <- factor(gene_ordering_file$V2, levels=c(paste0("chr", seq(1:23)), "chrX", "chrY", "chrM")) 
gene_ordering_file <- gene_ordering_file[order(gene_ordering_file$V2), ]
gene_ordering_file$V2 <- as.character(gene_ordering_file$V2)                 
write.table(gene_ordering_file,"/path/to/InferCNV/gene_ordering_file.txt"), sep="\t", col.names = F, row.names = F, quote=F)


### Create infercnv dataset for steroid cells

Idents(dataset) <- dataset$celltype_clusters
dataset_steroid <- subset(dataset, idents="Steroid cells")
Idents(dataset_steroid) <- dataset_steroid$orig.ident

# remove sample with < 200 cells
dataset_steroid <- subset(dataset_steroid, idents="ACC3", invert=T) 

counts_matrix <- dataset_steroid@assays$RNA@counts[ ,colnames(dataset_steroid)]
cell_annot <- cbind(names(dataset_steroid@active.ident), as.character(dataset_steroid@active.ident))
write.table(cell_annot,"/path/to/InferCNV/cell_annotations.txt"), sep="\t", col.names = F, row.names = F, quote=F)

infercnv_obj <- CreateInfercnvObject(raw_counts_matrix = counts_matrix,
                                     annotations_file = "/path/to/InferCNV/cell_annotations.txt",
                                     gene_order_file = "/path/to/InferCNV/gene_ordering_file.txt"),
                                     delim = "\t",
                                     ref_group_names = c("NAd1","NAd2","NAd3","NAd4"))
                    

### Run InferCNV

infercnv_obj <- infercnv::run(infercnv_obj,
                             cutoff = 0.1,  # use 0.1 for 10x-genomics
                             out_dir = "/path/to/InferCNV/outputs", 
                             cluster_by_groups = T,
                             analysis_mode = "samples",
                             denoise = T, 
                             HMM = F, 
                             num_threads = 12, 
                             useRaster = F)


### Classify cells using InferCNV

inferCNV_result <- readRDS("/path/to/InferCNV/outputs/run.final.infercnv_obj")
ref <- inferCNV_result@expr.data[, unlist(inferCNV_result@reference_grouped_cell_indices), drop=FALSE]
ref <- ref[, intersect(colnames(ref), rownames(dataset@meta.data))]
tum <- inferCNV_result@expr.data[, unlist(inferCNV_result@observation_grouped_cell_indices), drop=FALSE]
tum <- tum[, intersect(colnames(tum), rownames(dataset@meta.data))]

# Binary transformation for each  gene : normal(0)/altered(1) according to ref population distribution
CNA_value <- cbind(ref, tum)
for (i in 1:nrow(CNA_value)){
  CNA_value_gene <- CNA_value[i,]
  quantile_ref <- quantile(ref[i,], probs=seq(0,1,0.005))
  cut_low <- as.numeric(quantile_ref["0.5%"])
  cut_high <- as.numeric(quantile_ref["99.5%"])
  CNA_value_gene[which(CNA_value[i,] >= cut_low & CNA_value[i,] <= cut_high)] <- 0
  CNA_value_gene[which(CNA_value[i,] < cut_low | CNA_value[i,] > cut_high)] <- 1
  CNA_value[i,] <- CNA_value_gene }

# Compute proportion of genome altered
CNA_metrics <- as.data.frame(colSums(CNA_value)/nrow(CNA_value))

# Add annotations
annot <- cbind(as.character(dataset@meta.data$histotype), as.character(dataset@active.ident))
rownames(annot) <- rownames(dataset@meta.data)
Id_to_merge <- intersect(rownames(CNA_metrics), rownames(annot))
CNA_metrics <- cbind(CNA_metrics[Id_to_merge,], annot[Id_to_merge,])
colnames(CNA_metrics)[3:4]<- c("histotype", "celltype")
CNA_metrics$sample <- gsub("_.*" , "", rownames(CNA_metrics))

# Examine distribution of proportion of genome altered in cancer cells vs normal adrenals
boxplot(CNA_metrics$CNA_prop ~ CNA_metrics$histotype, outline=F)
abline(h=0.03, col="red")

# Classify cells based on proportion of genome altered
CNA_metrics$CNA_status <- rep(NA, nrow(CNA_metrics))
CNA_metrics$CNA_status[which(CNA_metrics$CNA_prop > 0.03)] <- "Tumor cells"
CNA_metrics$CNA_status[which(CNA_metrics$CNA_prop <= 0.03)] <- "Normal cells"

dataset$CNA_score <- rep(NA, nrow(dataset@meta.data))
dataset$CNA_score[rownames(CNA_metrics)] <-  CNA_metrics_all$CNA_prop
dataset$CNA_status <- rep(NA, nrow(dataset@meta.data))
dataset$CNA_status[rownames(CNA_metrics)]<-  CNA_metrics_all$CNA_status

saveRDS(dataset, "/path/to/Seurat/objects/Merge/dataset_filtered_normalized_annotated.RDS")


