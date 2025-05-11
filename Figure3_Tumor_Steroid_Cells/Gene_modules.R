library(ggplot2)
library(corrplot)
library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(tibble)
library(GSVA)
library(xlsx)
library(cowplot)
library(cola)
library(ggcorrplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

dataset <- readRDS("/path/to/Seurat_objects/Merge/dataset_filtered_normalized_annotated.RDS")


### Subset tumor steroid cells
# Exclude tumors with < 250 steroid cells

dataset_tumor <- dataset
DefaultAssay(dataset_tumor) <- "RNA"
dataset_tumor[['SCT']] <- NULL
Idents(dataset) <- dataset_tumor$orig.ident
dataset <- subset(dataset_tumor, idents = c("NAd1", "NAd2", "NAd3", "NAd4", "ACC3", "ACC14"), invert = T)
Idents(dataset_tumor) <- dataset_tumor$celltype_clusters
dataset_tumor <- subset(dataset_tumor, idents = "Steroid cells")
all_genes <- rownames(dataset_tumor)
dataset_list <-  SplitObject(dataset_tumor, split.by = "orig.ident")


### Run PCA for each patient

dataset_list <- lapply(X = dataset_list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = "percent.mt")
  x <- RunPCA(x, npcs = 20) })


### Select top PC genes (positive and negative) for each patient and merge

nPC = 10
n_top_genes = 50

list_genes <- list()
for (i in 1:32)
{
  dataset_tmp <- dataset_list[[i]]
  ech <- as.character(unique(dataset_tmp$orig.ident))
  mat_PCA <- dataset_tmp@reductions$pca@feature.loadings[,1:nPC]
  for (j in 1:nPC)
  {
    genes_pos <- names(sort(mat_PCA[,j], decreasing = T))[1:n_top_genes]
    genes_neg <- names(sort(mat_PCA[,j]))[1:n_top_genes]
    list_genes[[paste0(ech, "_PC", j, "_pos")]] <- genes_pos
    list_genes[[paste0(ech, "_PC", j, "_neg")]] <- genes_neg
  }
}


### Sorensen matrix 

sor_matrix <- matrix(nrow = length(list_genes), ncol = length(list_genes))
rownames(sor_matrix) <- names(list_genes)
colnames(sor_matrix) <- names(list_genes)
for (a in seq(ncol(sor_matrix))) {
  for (b in seq(ncol(sor_matrix))) {
    sor_matrix[a,b] = length(intersect(list_genes[[a]], list_genes[[b]])) / n_top_genes  } }

### Annotations for heatmaps

list_C1A <- c("ACC1a", "ACC1b", "ACC2a", "ACC2b", "ACC4", "ACC5", "ACC6", "ACC7", "ACC8a", "ACC8b")
list_C1B <- c("ACC9" ,"ACC10", "ACC11", "ACC12", "ACC13", "ACC15", "ACC16", "ACC17")

metadata <- data.frame(as.factor(gsub("_.*", "", names(list_genes))), row.names = names(list_genes))
colnames(metadata) <- "sample"
metadata$histotype <- rep(NA, nrow(metadata))
metadata$histotype[grep("ACA",metadata$sample)] <- "ACA"
metadata$histotype[grep("PBMAH",metadata$sample)] <- "PBMAH"
metadata$histotype[which(metadata$sample %in% list_C1A)] <- "ACC C1A"
metadata$histotype[which(metadata$sample %in% list_C1B)] <- "ACC C1B"
metadata$histotype <- as.factor(metadata$transcriptome)
list_color <- list(histotype = c("ACA" = "deepskyblue1", "PBMAH" = "mediumaquamarine", "ACC C1A" = "darkmagenta", "ACC C1B" = "hotpink3"))
col_sample <- c("#0D0887FF", "#20068FFF", "#2E0595FF", "#39049AFF","#5002A2FF", "#5B01A5FF", "#6600A7FF", "#7100A8FF", "#7B02A8FF", "#8606A6FF", "#8F0DA4FF", "#99159FFF", "#A21C9AFF", "#AB2394FF",
"#B22B8FFF", "#C13B82FF", "#C8437BFF", "#CE4B75FF", "#D5536FFF", "#DA5B69FF", "#E06363FF", "#E56B5DFF", "#E97357FF", "#EE7B51FF", "#F2844BFF", "#F68D45FF", "#F9973FFF", "#FBA139FF",
"#FDAB33FF", "#FDB52EFF", "#FEBF29FF", "#FDCA26FF")
names(col_sample) <- c(list_C1A, list_C1B, paste0("ACA", seq(1,8)), paste0("PBMAH", seq(1,6)))
list_color[["sample"]] <- col_sample


### Filter on Sorensen indexes

test <- apply(sor_matrix, 1, function(x) length(which(x >= 0.4)))
sor_matrix <- sor_matrix[which(test >= 3), which(test >= 3)]


### Clustering & defining gene modules
              
clust_sor <- pheatmap(sor_matrix, clustering_method = "ward.D2", 
              annotation_col = metadata,  annotation_colors = list_color, 
              show_rownames = F, show_colnames = F) 

ngroups = 8
groups <- cutree(a$tree_col, ngroups)
metadata[, "Gene modules"] <- rep(NA, nrow(metadata))
metadata[names(groups),"Gene modules"] <- groups
metadata[,"Gene modules"] <- factor(as.character(metadata[,"Gene modules"]), 
                                                 levels = c("1","2","3","4","5","6","7","8"), 
                                                 labels = c("GM1", "GM5", "GM2", "GM4", "GM3", "GM6", "GM7", "GM8"))

list_color[["Gene modules"]] <-  c("GM1" = "#7A0403FF", "GM2" = "#DB3A07FF" , "GM3" = "#FE9B2DFF" , "GM4" = "#62FC6BFF",
                                   "GM5" = "#D2E935FF", "GM6" = "#1BD0D5FF", "GM7" = "#30123BFF",  "GM8" = "#4777EFFF")

list_prog <- NULL
for (group in 1 : length(unique(groups)))
{
  sel <- names(groups)[which(groups == group)]
  list_prog_group <- list_genes[sel]
  samples <- unique(gsub("_.*", "", sel))
  
  # unique genes by sample
  list_prog_by_sample <- list() 
  for (sample in samples)
  {
    genes_sample <- unique(unlist(list_prog_group[grep(sample, names(list_prog_group))]))
    list_prog_by_sample[[sample]] <- genes_sample
  }
  
  # keep genes present in programs from >= 50% ech or >=2 samples
  table_recurrence <- sort(table(unlist(list_prog_by_sample)), decreasing=T)
  n_sample <- length(list_prog_by_sample)
  genes_sel <- names(table_recurrence)[which(table_recurrence >= max(2, n_ech/2))]
  list_prog[[paste0("prog",group)]] <- genes_sel
}

names(list_prog) <- c("GM1_ZG", "GM5_ECM", "GM2_ZF1", "GM4_ZR", "GM3_ZF2", "GM6_Translation", "GM7_Mitosis", "GM8_Hypoxia")


### Overrepresentation analysis

# define genesets and universes
df_go <- as.data.frame(org.Hs.egGO)
go_universe <- unique(sort(df_go$gene_id))

# create geneset normal adrenal
Nad_geneset <- read.csv("/path/to/SupplTable3.csv")
clusters <- as.character(unique(Nad_geneset$cluster))
EntrezID <- bitr(Nad_geneset$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
EntrezID <- EntrezID[!duplicated(EntrezID$SYMBOL),]
Nad_geneset$ENTREZID <- rep(NA, nrow(Nad_geneset))
for (ID in EntrezID$SYMBOL) {
  Nad_geneset$ENTREZID[which(Nad_geneset$gene == ID)] <- EntrezID[which(EntrezID$SYMBOL == ID), "ENTREZID"] }
Nad_geneset <- Nad_geneset[which(!is.na(Nad_geneset$ENTREZID)),]
gs2gene_adrenal <- as.data.frame(cbind( gsub(" ", "_", Nad_geneset$cluster), Nad_geneset$ENTREZID))
colnames(gs2gene_adrenal) <- c("gs_name", "entrez_gene")
gs2gene_adrenal <- tibble(gs2gene_adrenal)
gs2name_adrenal <- as.data.frame( cbind( gsub(" ", "_", Nad_geneset$cluster), Nad_geneset$gene))
colnames(gs2name_adrenal) <-  c("gs_name", "gene_symbol")
gs2name_adrenal <- tibble(gs2name_adrenal)
gs2gene_adrenalgo <- rbind(as.data.frame(gs2gene_adrenal), as.data.frame(gs2gene_Go))
gs2gene_adrenalgo <- tibble(gs2gene_adrenalgo)
gs2name_adrenalgo <- rbind(as.data.frame(gs2name_adrenal), as.data.frame(gs2name_Go))
gs2name_adrenalgo <- tibble(gs2name_adrenalgo)

# convert gene ID for each cluster
gene_modules <- names(list_prog)
GM_markers <- NULL
for (gm in gene_modules)
{
  genes_gm <- list_prog[[gm]]
  name_cluster <- gm
  genes_EntrezID <- bitr(genes_gm, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
  gm_markers <- data.frame(EntrezID = genes_EntrezID$ENTREZID, clusterID = name_cluster)
  GM_markers <- rbind(GM_markers, gm_markers) }
}

# compare clusters
Nad_enricher <- compareCluster(EntrezID ~ clusterID,  data = metaprog_markers, fun = "enricher", TERM2GENE = gs2gene_adrenal, TERM2NAME = gs2name_adrenal)
Nad_enricher@compareClusterResult$Description <- Nad_enricher@compareClusterResult$ID
Go_enricher <-  compareCluster(EntrezID ~ clusterID,  data = metaprog_markers, fun = "enricher", TERM2GENE = gs2gene_Go, TERM2NAME = gs2name_Go)
Go_enricher@compareClusterResult$Description <- Go_enricher@compareClusterResult$ID
combined_enricher <- Nad_enricher
combined_enricher@compareClusterResult <- rbind(combined_enricher@compareClusterResult, Go_enricher@compareClusterResult) 
combined_enricher@compareClusterResult <- combined_enricher@compareClusterResult[order(combined_enricher@compareClusterResult$clusterID), ]
             
              
### Figure 3

# Heatmap
p_heatmap <- pheatmap(sor_matrix, clustering_method = "ward.D2", 
                      annotation_col = metadata,  annotation_colors = list_color, 
                      show_rownames = F, show_colnames = F)

# DotPlot Enrichment analyses
dotplot(combined_enricher, showCategory = 3, font.size = 8)




