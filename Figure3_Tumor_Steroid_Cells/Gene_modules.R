library(ggplot2)
library(corrplot)
library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(tibble)
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
for (i in 1:32) {
  dataset_tmp <- dataset_list[[i]]
  ech <- as.character(unique(dataset_tmp$orig.ident))
  mat_PCA <- dataset_tmp@reductions$pca@feature.loadings[,1:nPC]
  for (j in 1:nPC) {
    genes_pos <- names(sort(mat_PCA[,j], decreasing = T))[1:n_top_genes]
    genes_neg <- names(sort(mat_PCA[,j]))[1:n_top_genes]
    list_genes[[paste0(ech, "_PC", j, "_pos")]] <- genes_pos
    list_genes[[paste0(ech, "_PC", j, "_neg")]] <- genes_neg } }


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

list_GM <- NULL
for (group in 1 : length(unique(groups))) {
  sel <- names(groups)[which(groups == group)]
  list_GM_group <- list_genes[sel]
  samples <- unique(gsub("_.*", "", sel))
  # unique genes by sample
  list_GM_by_sample <- list() 
  for (sample in samples) {
    genes_sample <- unique(unlist(list_GM_group[grep(sample, names(list_GM_group))]))
    list_GM_by_sample[[sample]] <- genes_sample }
  # keep genes present in programs from >= 50% ech or >=2 samples
  table_recurrence <- sort(table(unlist(list_GM_by_sample)), decreasing=T)
  n_sample <- length(list_GM_by_sample)
  genes_sel <- names(table_recurrence)[which(table_recurrence >= max(2, n_ech/2))]
  list_GM[[paste0("prog",group)]] <- genes_sel }

names(list_GM) <- c("GM1_ZG", "GM5_ECM", "GM2_ZF1", "GM4_ZR", "GM3_ZF2", "GM6_Translation", "GM7_Mitosis", "GM8_Hypoxia")
saveRDS(list_GM, "/path/to/GeneModules/results/Gene_modules.RDS")


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
             

### Add gene modules scores to Sn-RNAseq dataset

# Add module score
for (i in 1:length(list_prog)) {
  program <- list_prog[i]
  dataset <- AddModuleScore(dataset, features = program, name = names(program)) }
col_to_replace <- grep("GM", colnames(dataset@meta.data))[1:8]
for (i in col_to_replace) { colnames(dataset@meta.data)[i] <-  substr(colnames(dataset@meta.data)[i], 1, (nchar(colnames(dataset@meta.data)[i])-1) ) }

# Assign cells to metaprogramme max
test <- apply(dataset@meta.data[,names(list_prog)], 1, function(x) which(x == max(x)))
dataset$GM_maj <- names(list_prog)[test]
dataset$GM_maj[which(dataset$celltype_clusters != "Steroid cells")] <- NA

saveRDS("/path/to/Seurat_objects/Merge/dataset_filtered_normalized_annotated_GM.RDS")
              
### Figure 3 plots

# Heatmap
p_heatmap <- pheatmap(sor_matrix, clustering_method = "ward.D2", 
                      annotation_col = metadata,  annotation_colors = list_color, 
                      show_rownames = F, show_colnames = F)

# DotPlot Enrichment analyses
dotplot(combined_enricher, showCategory = 3, font.size = 8)

# UMAP by gene module
GM_colors <-  c( "#7A0403FF", "#DB3A07FF", "#FE9B2DFF", "#62FC6BFF", "#D2E935FF", "#1BD0D5FF", "#30123BFF", "#4777EFFF")
p_umap_GM <- DimPlot(dataset, reduction = "umap", group.by = "GM_maj", raster = F, 
        cols=GM_colors, na.value = "grey70") +
  ggtitle("") 
              
# proportion of cells in each sample 
dataset$orig.ident <- factor(dataset$orig.ident, 
                             levels = c("ACC1a","ACC1b","ACC2a","ACC2b","ACC3","ACC4","ACC5","ACC6","ACC7","ACC8a","ACC8b",
                                        "ACC9","ACC10","ACC11","ACC12","ACC13","ACC14","ACC15","ACC16","ACC17",
                                        "ACA1","ACA2","ACA3","ACA4","ACA5","ACA6","ACA7","ACA8",
                                        "PBMAH1","PBMAH2","PBMAH3","PBMAH4","PBMAH5","PBMAH6"))
Idents(dataset) <- dataset$celltype_clusters
dataset_test <- subset(dataset, idents = "Steroid cells")
Idents(dataset_test) <- dataset_test$celltype_ScibetGarnett
dataset_test <- subset(dataset_test, idents = "Steroid cells")
temp_Names<- unique(as.character(dataset_test$orig.ident))
temp_colname<- "GM_maj"
temp_colour_pal <- data.frame(metaprog_maj = c("GM1_ZG", "GM2_ZF1", "GM3_ZF2", "GM4_ZR", "GM5_ECM", "GM6_Translation", "GM7_Mitosis", "GM8_Hypoxia"),
                                colour = GM_colors))
temp_cellprop_df_all <- NULL
for(i in c(1:length(temp_Names))) 
{
  # i=1
  temp_filtered_summary <- subset(dataset_test@meta.data, orig.ident == temp_Names[i])
  temp_cellprop_df <- data.frame(unclass(table(temp_filtered_summary[,temp_colname])))
  temp_cellprop_df_all <- rbind(temp_cellprop_df_all,
                                (temp_cellprop_df <- data.frame(sample = rep(temp_Names[i], times = length(row.names(temp_cellprop_df))),
                                                                metaprog_maj = row.names(temp_cellprop_df), value = temp_cellprop_df[,1],
                                                                proportions = temp_cellprop_df[,1]/colSums(temp_cellprop_df),
                                                                subtype = unique(temp_filtered_summary$histotype))))
}
temp_cellprop_df_all$metaprog_maj <- factor(temp_cellprop_df_all$GM_maj, levels=temp_colour_pal$GM_maj)
temp_cellprop_df_all$subtype <- factor(temp_cellprop_df_all$subtype, levels=c("ACC C1A","ACC C1B","ACA", "PBMAH", "Normal adrenal"))
p_prop <- ggplot(temp_cellprop_df_all, aes(x=sample, fill=GM_maj, y=proportions)) +
  geom_bar(stat="identity", color = "black") +
  theme(axis.text.x=element_text(angle=45, size=7, hjust=1),
        strip.background = element_rect(colour="black", fill="white", linewidth=0.5, linetype="solid"),
        strip.text.x = element_text(size=8, face = "bold"),
        axis.text = element_text(size=7), 
        axis.title = element_text(size=8, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.key.size = unit(0.5, 'cm'), 
        legend.key.height = unit(0.5, 'cm'), 
        legend.key.width = unit(0.5, 'cm'), 
        legend.title = element_text(size=7),
        legend.text = element_text(size=7)) + 
  xlab("Patient ID") +
  ylab("GM Proportions among Steroid cells")  +
  guides(fill=guide_legend(title="Gene modules")) +
  scale_fill_manual(values = as.vector(temp_colour_pal$colour)) +
  facet_grid(. ~ subtype, scales="free", space = "free")





