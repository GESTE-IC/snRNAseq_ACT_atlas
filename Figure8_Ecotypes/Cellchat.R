library(Matrix)
library(dplyr)
library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(tidyverse)
library(viridis)
library(CellChat)
library(patchwork)
library(pheatmap)

dataset <- readRDS("/path/to/Seurat_objects/Merge/dataset_filtered_normalized_annotated_GM_TME.RDS")
ecotypes <- readRDS("path/to/Ecotypes/results/Ecotypes_consensus.RDS")


### prepare annotations and datasets for cell-cell interactions analyses

# Add CCI annotations
dataset$celltype_CCI <- as.character(dataset$celltype_deconv)
dataset$celltype_CCI[which(dataset$celltype_CCI == "Steroid_cells")] <- dataset[which(dataset$celltype_CCI == "Steroid_cells"), "GM_maj"]
Idents(dataset) <- dataset$celltype_CCI
dataset <-  subset(dataset, idents = names(ecotypes))

# list all possible CCI within ecotypes
possible_CCI <- NULL
table_CCs <- data.frame(pop_1_2=NULL, ecotype=NULL)
for (eco in unique(ecotypes)) {
  celltypes <- names(ecotypes)[which(ecotypes == eco)]
  list_CC <- NULL
  for (celltype in celltypes) { list_CC <- c(list_CC, paste(celltype, celltypes, sep="_")) }
  pop_1_2 <- unique(list_CC)
  ecotype <- rep(eco, length(list_CC))
  table_CCs <- rbind(table_CCs, cbind(pop_1_2, ecotype))
  possible_CCI <- c(possible_CCI, list_CC) }
rownames(table_CCs) <- table_CCs$pop_1_2


### Data input & processing and initialization of CellChat object

# Prepare dataset
DefaultAssay(dataset) <- "RNA"
cellchat <- createCellChat(object = dataset, meta = dataset@meta.data, group.by = "celltype_CCI")
cellchat <- setIdent(cellchat, ident.use = "celltype_CCI")

# Prepare database
CellChatDB <- CellChatDB.human 
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

# Preprocess expression data
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)


### Inference of cell-cell communication network

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "triMean", trim= 0.1, raw.use = FALSE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

# Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Calculate and vizualize the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


### Table of interactions

df_all_prob <- data.frame(pathway=NA, interaction=NA, ligand=NA, receptor=NA, pop_cell=NA, prob=NA, pval=NA)
for (i in 1:dim(cellchat@net$prob)[3]) {
  prob_res <- cellchat@net$prob[,,i]
  pval_res <- cellchat@net$pval[,,i]
  for (a in 1:nrow(prob_res)) {
    cell_pop <- paste(rownames(prob_res)[a], colnames(prob_res), sep="_")
    interaction <- rep(cellchat@LR$LRsig$interaction_name_2[i], nrow(prob_res))
    pathway <- rep(cellchat@LR$LRsig$pathway_name[i], nrow(prob_res))
    ligand <- rep(cellchat@LR$LRsig$ligand[i], nrow(prob_res))
    receptor <- rep(cellchat@LR$LRsig$receptor[i], nrow(prob_res))
    prob <- prob_res[a,]
    pval <- pval_res[a,]
    df_prob <- data.frame(pathway, interaction, ligand, receptor, pop_cell, prob, pval)
    df_all_prob <- rbind(df_all_prob, df_prob) }
}
saveRDS(df_all_prob, "/path/to/Ecotypes/results/Table_interactions_Cellchat.RDS"))

# filtering signif interactions
df_all_prob <- df_all_prob[which(df_all_prob$pval<0.05),]

# filtering by ecotype
df_all_prob <- df_all_prob[which(df_all_prob$cell_pop %in% possible_CCI), ]
df_all_prob$ecotype <- rep(NA, nrow(df_all_prob))
for (i in 1:nrow(df_all_prob)) {
  popCCI <- df_all_prob[i,"cell_pop"]
  df_all_prob[i,"ecotype"] <- table_CCs[which(table_CCs$pop_1_2 == popCCI), "ecotype"]
}

# ordering by pathway, by ecotype and by pval
df_all_prob <- df_all_prob[order(df_all_prob$pval),]
df_all_prob <- df_all_prob[order(df_all_prob$ecotype),]
df_all_prob <- df_all_prob[order(df_all_prob$prob, decreasing=T),]
df_all_prob <- df_all_prob[order(df_all_prob$pathway),]


### Table of pathways

# create table
df_all_prob_pathway <- data.frame(pathway=NA, pop_1_2=NA, pop_1_type=NA, pop_2_type=NA, prob=NA)
for (i in 1:length(pathways.show.all)) {
  pathway_res <- cellchat@netP$prob[,,i]
  for (a in 1:nrow(pathway_res)) {
    pop_1_2 <- paste(rownames(pathway_res)[a], colnames(pathway_res), sep="_")
    pop_1_type <- rownames(pathway_res)[a]
    pop_2_type <- colnames(pathway_res)
    pathway <- rep(pathways.show.all[i], nrow(pathway_res))
    prob <- pathway_res[a,]
    df_prob <- data.frame(pathway, pop_1_2, pop_1_type, pop_2_type, prob)
    df_all_prob_pathway <- rbind(df_all_prob_pathway, df_prob) }
}
df_all_prob_pathway <- df_all_prob_pathway[which(!is.na(df_all_prob_pathway$prob)), ]
df_all_prob_pathway <- df_all_prob_pathway[order(df_all_prob_pathway$prob, decreasing = T),]
df_all_prob_pathway <- df_all_prob_pathway[order(df_all_prob_pathway$pathway),]

# add cell type
df_all_prob_pathway$pop_1_type[grep("GM", df_all_prob_pathway$pop_1_type)] <- "Steroid"
df_all_prob_pathway$pop_1_type[which(df_all_prob_pathway$pop_1_type %in% c("Resident_macrophages_1", "TAM1", "Perivascular_TAM", "Inflammatory_macrophages"))] <- "Myelo"
df_all_prob_pathway$pop_1_type[which(df_all_prob_pathway$pop_1_type %in% c( "Naive_Memory_T_cells", "Exhausted_T_cells", "NK_like_T_cells"))] <- "Lympho"
df_all_prob_pathway$pop_1_type[which(df_all_prob_pathway$pop_1_type %in% c( "CAF1", "Resident_fibroblasts_1", "CAF3"))] <- "Fibro"
df_all_prob_pathway$pop_1_type[which(df_all_prob_pathway$pop_1_type %in% c( "EC_venous", "EC_arterial", "TEC2"))] <- "Endoth"
df_all_prob_pathway$pop_2_type[grep("GM", df_all_prob_pathway$pop_2_type)] <- "Steroid"
df_all_prob_pathway$pop_2_type[which(df_all_prob_pathway$pop_2_type %in% c("Resident_macrophages_1", "TAM1", "Perivascular_TAM", "Inflammatory_macrophages"))] <- "Myelo"
df_all_prob_pathway$pop_2_type[which(df_all_prob_pathway$pop_2_type %in% c( "Naive_Memory_T_cells", "Exhausted_T_cells", "NK_like_T_cells"))] <- "Lympho"
df_all_prob_pathway$pop_2_type[which(df_all_prob_pathway$pop_2_type %in% c( "CAF1", "Resident_fibroblasts_1", "CAF3"))] <- "Fibro"
df_all_prob_pathway$pop_2_type[which(df_all_prob_pathway$pop_2_type %in% c( "EC_venous", "EC_arterial", "TEC2"))] <- "Endoth"
df_all_prob_pathway$pop_1_2_type <- paste(df_all_prob_pathway$pop_1_type, df_all_prob_pathway$pop_2_type, sep="_")

# scale by cell pop
df_all_prob_scale <- matrix(nrow=0, ncol=ncol(df_all_prob_pathway))
colnames(df_all_prob_scale) <- colnames(df_all_prob_pathway)
all_couples <- unique(df_all_prob_pathway$pop_1_2_type)
all_pathways <- unique(df_all_prob_pathway$pathway)
for (couple in all_couples) {
  df_couple <- df_all_prob_pathway[which(df_all_prob_pathway$pop_1_2_type == couple ), ]
  df_couple_new <- matrix(nrow=0, ncol=ncol(df_couple))
  colnames(df_couple_new) <- colnames(df_couple)
  for (pathway in all_pathways) {
    # pathway = "VEGF"
    df_couple_pathway <- df_couple[which(df_couple$pathway == pathway), ]
    if (sd(df_couple_pathway$prob) > 0 ) {
    df_couple_pathway$prob <- scale(df_couple_pathway$prob) 
    df_couple_new <- rbind(df_couple_new, df_couple_pathway) }
    if (sd(df_couple_pathway$prob) == 0 ) {
      df_couple_new <- rbind(df_couple_new, df_couple_pathway) }  }
    df_all_prob_scale <- rbind(df_all_prob_scale, df_couple_new) }

# Table of results fitered and annotated 
df_all_prob_ordered <- df_all_prob_pathway[rownames(df_all_prob_scale), ]
df_all_prob_filtered <- df_all_prob_ordered[which(df_all_prob_scale$prob >2), ]
# df_all_prob_filtered <- df_all_prob_ordered[which(df_all_prob_scale$prob >3), ]
df_all_prob_filtered <- df_all_prob_filtered[which(df_all_prob_filtered$pop_1_2 %in% possible_CCI), ]
df_all_prob_filtered$ecotype_enriched <- rep(NA, nrow(df_all_prob_filtered))
df_all_prob_filtered[which(df_all_prob_filtered$pop_1_2 %in% table_CCs[which(table_CCs$ecotype == "Eco1"),"pop_1_2"]), "ecotype_enriched"] <- "Eco1"
df_all_prob_filtered[which(df_all_prob_filtered$pop_1_2 %in% table_CCs[which(table_CCs$ecotype == "Eco2"),"pop_1_2"]), "ecotype_enriched"] <- "Eco2"
df_all_prob_filtered[which(df_all_prob_filtered$pop_1_2 %in% table_CCs[which(table_CCs$ecotype == "Eco3"),"pop_1_2"]), "ecotype_enriched"] <- "Eco3"
df_all_prob_filtered[which(df_all_prob_filtered$pop_1_2 %in% table_CCs[which(table_CCs$ecotype == "Eco4"),"pop_1_2"]), "ecotype_enriched"] <- "Eco4"
df_all_prob_filtered[which(df_all_prob_filtered$pop_1_2 %in% table_CCs[which(table_CCs$ecotype == "Eco5"),"pop_1_2"]), "ecotype_enriched"] <- "Eco5"
df_all_prob_filtered$pop_1_type <- NULL
df_all_prob_filtered$pop_2_type <- NULL
colnames(df_all_prob_filtered) <- c("pathway", "cellstate1_cellstate2", "prob", "cellpop1_cellpop2", "ecotype_enriched")
df_all_prob_filtered <- df_all_prob_filtered[,  c("pathway", "cellstate1_cellstate2", "cellpop1_cellpop2", "prob", "ecotype_enriched")]
df_all_prob_filtered <- df_all_prob_filtered[order(df_all_prob_filtered$pathway), ]
df_all_prob_filtered <- df_all_prob_filtered[order(df_all_prob_filtered$ecotype_enriched), ]


### Heatmap by pathway and ecotype

# create dataframe of interaction scaled by pathway
df_by_pathway <- matrix(nrow=length(unique(df_all_prob_scale$pop_1_2)), ncol=length(unique(df_all_prob_scale$pathway)))
colnames(df_by_pathway) <- unique(df_all_prob_scale$pathway)
rownames(df_by_pathway) <- unique(df_all_prob_scale$pop_1_2)
for (pathway in unique(df_all_prob_scale$pathway))   {
    names_pop_cell_ordered <- df_all_prob_scale[which(df_all_prob_scale == pathway), "pop_1_2"]
    res_pop_cell_ordered <- df_all_prob_scale[which(df_all_prob_scale == pathway), "prob"]
    df_by_pathway[names_pop_cell_ordered, pathway] <- res_pop_cell_ordered }

# filter on sd
df_by_pathway_sel <- df_by_pathway[intersect(possible_CCI, rownames(df_by_pathway)), ]
test1 <- apply(df_by_pathway_sel, 1, function(x) length(which(x>3)))
df_by_pathway_sel <- df_by_pathway_sel[names(test1)[which(test1>0)], ]
test2 <- apply(df_by_pathway_sel, 2, function(x) length(which(x>3)))
df_by_pathway_sel <- df_by_pathway_sel[, names(test2)[which(test2>0)]]

# annotations for heatmap
annot_eco <- rep(NA, nrow(df_by_pathway_sel))
names(annot_eco) <- rownames(df_by_pathway_sel)
annot_eco[which(names(annot_eco) %in% table_CCs[which(table_CCs$ecotype == "Eco1"), "pop_1_2"])] <- "Eco1"
annot_eco[which(names(annot_eco)  %in% table_CCs[which(table_CCs$ecotype == "Eco2"), "pop_1_2"])] <- "Eco2"
annot_eco[which(names(annot_eco)  %in% table_CCs[which(table_CCs$ecotype == "Eco3"), "pop_1_2"])] <- "Eco3"
annot_eco[which(names(annot_eco)  %in% table_CCs[which(table_CCs$ecotype == "Eco4"), "pop_1_2"])] <- "Eco4"
annot_eco[which(names(annot_eco)  %in% table_CCs[which(table_CCs$ecotype == "Eco5"), "pop_1_2"])] <- "Eco5"
annot <- data.frame(annot_eco)
colnames(annot) <- "Ecotype"
list_color <- list(Ecotype = c("Eco1"="#3f2d54", "Eco2"="#d18975", "Eco3"="#2d543d", "Eco4"="#8fd175", "Eco5"="#75b8d1"))

# color scaling
df_by_pathway_sel[which(df_by_pathway_sel>4)] <- 4

# figure           
p_ecotypes_CCI <- pheatmap(df_by_pathway_sel,  clustering_method = "ward.D2",
           annotation_row = annot, annotation_colors = list_color, angle_col = 45)

# Suppl Table 18
df_all_prob <- df_all_prob[which(df_all_prob$pathway %in% colnames(df_by_pathway_sel)), ]


### Visualization of CCI in selected pathways

pathways_Eco5 <- c("FASLG", "LT", "IL1","IL4", "IL6", "MSTN", "TGFb", "COMPLEMENT", "CCL", "CXCL", "CD40", "CD80", "CD86", "ICAM")
pathways_Eco2 <- c("VCAM", "ALCAM", "CD6", "CDH", "CD99", "JAM")
pathways_Eco1 <- c("VEGF", "NECTIN", "TENASCIN", "EPO", "HGF")

p_Eco1 <- netVisual_chord_gene(cellchat, sources.use = c(2,12,13,21), targets.use = c(2,12,13,21), signaling = pathways_Eco1, lab.cex = 1.3, legend.pos.y = 30, reduce = 0.002)
p_Eco2 <- netVisual_chord_gene(cellchat, sources.use = c(5,6,8,10,11,18), targets.use = c(5,6,8,10,11,18), signaling = pathways_Eco2, lab.cex = 1.3, legend.pos.y = 30, reduce = 0.001)
p_Eco5 <- netVisual_chord_gene(cellchat, sources.use = c(1,9,14), targets.use = c(1,9,14), signaling = pathways_Eco5, lab.cex = 1.2, legend.pos.y = 30, reduce = 0.005)


