library(ggplot2)
library(corrplot)
library(Seurat)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(future)
library(tibble)
library(GSVA)
library(xlsx)
library(cowplot)
library(cola)
library(ggcorrplot)
library(survival)
library(survcomp)
library(ggpubr)
library(MASS)
library(lmtest)


### Load annotations and results of deconvolution

annot <- readRDS("/path/to/clinical/annotations/clinical_annotation_3cohorts.RDS")

res_GM <- readRDS("/path/to/Deconvolution/results/Deconvolution_GM.RDS")

res_cibersort_TCGA <-  read.csv2("/path/to/cibersort/results/TCGA/CIBERSORTx_Adjusted.txt", sep="\t")
rownames(res_cibersort_TCGA) <- res_cibersort_TCGA$Mixture
celltypes <- setdiff(colnames(res_cibersort_TCGA), c("Mixture","P.value","Correlation","RMSE"))
res_cibersort_TCGA[,celltypes] <- apply(res_cibersort_TCGA[,celltypes], 2, function(x) as.numeric(x))
                                       
res_cibersort_CIT <-  read.csv2("/path/to/cibersort/results/CIT/CIBERSORTx_Adjusted.txt", sep="\t")
rownames(res_cibersort_CIT) <- res_cibersort_CIT$Mixture
res_cibersort_CIT[,celltypes] <- apply(res_cibersort_CIT[,celltypes], 2, function(x) as.numeric(x))
                                        
res_cibersort_FFPE <-  read.csv2("/path/to/cibersort/results/FFPE/CIBERSORTx_Adjusted.txt", sep="\t")
rownames(res_cibersort_FFPE) <- res_cibersort_FFPE$Mixture
res_cibersort_FFPE[,celltypes] <- apply(res_cibersort_FFPE[,celltypes], 2, function(x) as.numeric(x))


 ### Aggregate and transform into Z scores                                       
                                        
res_TCGA <- cbind(res_cibersort_TCGA[ , celltypes], res_GM[rownames(res_cibersort_TCGA) , ])
res_Z_TCGA <- scale(res_TCGA)
                                        
res_CIT <- cbind(res_cibersort_CIT[ , celltypes], res_GM[rownames(res_cibersort_CIT) , ])
res_Z_CIT <- scale(res_CIT)
                                        
res_FFPE <- cbind(res_cibersort_FFPE[ , celltypes], res_GM[rownames(res_cibersort_FFPE) , ])
res_Z_FFPE <- scale(res_FFPE)

res_Z_all <- rbind(Z_CIT, Z_TCGA, Z_FFPE)
res_Z_all <- res_Z_all[, which(colnames(res_Z_all) != "Steroid_cells")]

                                        
### Association of single cell signatures with survival  -  Suppl Table 7  

# univariate analyses                                        
Varqual <- c("age", "sex", "cortisolsecretion", "stage", "mitoses", "ki67")
Varquant <- colnames(res_Z_all)
datatest <- cbind(annot[rownames(res_Z_all),], res_Z_all)                                     

ressurv <- matrix(nrow=0, ncol=8)
colnames(ressurv) <- c("variable", "N", "HR_DFS", "CI_DFS", "p_DFS", "HR_OS", "CI_OS", "p_OS")

for (v in Varqual)
{
  N <- length(which(!is.na(datatest[, v])))
  categ <- unique(datatest[, v])
  categ <- categ[!is.na(categ)]
  coxDFS <- summary(coxph(Surv(datatest$follow_up_DFS_months, datatest$DFS_status) ~ datatest[,v]))
  coxOS <- summary(coxph(Surv(datatest$follow_up_OS_months, datatest$OS_status) ~ datatest[,v]))
  tmp <- matrix(nrow = length(categ), ncol = 8)
  colnames(tmp) <- colnames(ressurv)
  tmp[1, "variable"] <- v
  tmp[2 : nrow(tmp), "variable"] <- paste("   ", categ[2:length(categ)], "vs", categ[1])
  tmp[1, "N"] <- N
  tmp[2:nrow(tmp), "HR_DFS"] <- round(coxDFS$coefficients[,2],2)
  tmp[2:nrow(tmp), "CI_DFS"] <- paste0(round(coxDFS$conf.int[,3],2),"-",round(coxDFS$conf.int[,4],2))
  tmp[2:nrow(tmp), "p_DFS"] <- signif(coxDFS$coefficients[,5],3)
  tmp[2:nrow(tmp), "HR_OS"] <- round(coxOS$coefficients[,2],2)
  tmp[2:nrow(tmp), "CI_OS"] <- paste0(round(coxOS$conf.int[,3],2),"-",round(coxOS$conf.int[,4],2))
  tmp[2:nrow(tmp), "p_OS"] <- signif(coxOS$coefficients[,5],3)
  ressurv <- rbind(ressurv, tmp)
}

for (v in Varquant)
{  
  N <- length(which(!is.na(datatest[, v])))
  coxDFS <- summary(coxph(Surv(datatest$follow_up_DFS_months, datatest$DFS_status) ~ datatest[ ,v]))
  coxOS <- summary(coxph(Surv(datatest$follow_up_OS_months, datatest$OS_status) ~ datatest [,v]))
  tmp <- matrix(nrow = 1, ncol = 8)
  colnames(tmp) <- colnames(ressurv)
  tmp[1, "N"] <- N
  tmp[, "variable"] <- paste(v, "Z-score (/1 unit increase))")
  tmp[, "HR_DFS"] <- round(coxDFS$coefficients[,2],2)
  tmp[, "CI_DFS"] <- paste0(round(coxDFS$conf.int[,3],2),"-",round(coxDFS$conf.int[,4],2))
  tmp[, "p_DFS"] <- signif(coxDFS$coefficients[,5],3)
  tmp[, "HR_OS"] <- round(coxOS$coefficients[,2],2)
  tmp[, "CI_OS"] <- paste0(round(coxOS$conf.int[,3],2),"-",round(coxOS$conf.int[,4],2))
  tmp[, "p_OS"] <- signif(coxOS$coefficients[,5],3)
  ressurv <- rbind(ressurv, tmp)
  
  datatest$median_score <- as.character(cut(datatest[,v], breaks=c((-1000), median(datatest[,v], na.rm=T), 1000), labels=c("low", "high")))
  pdf(paste0("/path/to/deconvolution/results/stats/KM_", v, "_DFS.pdf"), height = 7, width = 8)
  km.coxph.plot(Surv(follow_up_DFS_months,DFS_status)~median_score,data=datatest,x.label="Time from diagnosis (months)",y.label="Progression-free survival (%)",
                yscale=100,leg.text=paste(c("high", "low")," ",sep=""),leg.pos="bottomright",leg.bty='n', mark.time=T, 
                main.title="",show.n.risk=T, n.risk.step=12,xlim=c(0,120),.col=c('royalblue4','goldenrod1'),frame=F,.lty=1,.lwd=2,n.risk.cex=1)
  dev.off()
  pdf(paste0("/path/to/deconvolution/results/stats/KM_", v, "_OS.pdf"), height = 7, width = 8)
  km.coxph.plot(Surv(follow_up_OS_months,OS_status)~median_score,data=datatest,x.label="Time from diagnosis (months)",y.label="Overall survival (%)",
                yscale=100,leg.text=paste(c("high", "low")," ",sep=""),leg.pos="bottomright",leg.bty='n', mark.time=T, 
                main.title="",show.n.risk=T,n.risk.step=12,xlim=c(0,120),.col=c('royalblue4','goldenrod1'),frame=F,.lty=1,.lwd=2,n.risk.cex=1)
  dev.off()
}
                  

# multivariate analyses

complete_model_OS <- coxph(Surv(follow_up_OS_months, OS_status) ~ CAF1 + EC_arterial + Resident_fibroblasts_1 +
                          CAF3 + TEC2 + Exhausted_T_cells + Inflammatory_macrophages + 
                          GM3_ZF2 + GM7_Mitosis + stage, 
                        data = datatest)
step_model_OS <- stepAIC(complete_model_OS, direction = "both")                                        

complete_model_DFS <- coxph(Surv(follow_up_DFS_months, DFS_status) ~ CAF1 + EC_arterial + Resident_fibroblasts_1 +
                          CAF3 + TEC2 + Inflammatory_macrophages + 
                          GM3_ZF2 + GM7_Mitosis + stage, 
                        data = datatest)
step_model_DFS <- stepAIC(complete_model_DFS, direction = "both")                                 


###  Definition of ecotypes

# hierarchical clustering                                    
Test <- res_Z_all
Test[which(Test < (-3))] <- (-3)
Test[which(Test > 3)] <- 3
hclut_Test <- pheatmap(Test, clustering_method = "ward.D2", show_rownames = F, show_colnames = T, angle_col = 45) 

# consensus clustering annotations                                     
list_anno_col <- c("TAM1"="violetred4", "EC_venous"="lightcoral", "CAF1"="slateblue4", "CAF3"="slateblue3", "EC_arterial"="firebrick", "Resident_fibroblasts_1"="slateblue1",       
                   "Naive_Memory_T_cells"="darkgoldenrod", "Exhausted_T_cells"="darkgoldenrod4", "TEC2"="firebrick1", "Perivascular_TAM"="violetred1", "Inflammatory_macrophages"="violet",
                   "Resident_macrophages_1"="plum1", "NK_like_T_cells"="darkgoldenrod1", "GM1_ZG"="turquoise", "GM2_ZI_ZF"="turquoise4", "GM3_ZF1"="seagreen3", "GM4_ZR"="yellowgreen", "GM5_ZF2"="seagreen4",                 
                   "GM6_translation"="purple4", "GM7_Mitosis"="orchid4", "GM8_Hypoxia"="indianred")

# consensus clustering plot                                       
res <- consensus_partition(Test, max_k = 10, partition_repeat = 50, 
                           anno = colnames(Test), anno_col = list_anno_col, partition_method = "hclust",
                           sample_by="column")                                                                           
select_partition_number(res)
suggest_best_k(res)
p_consensus <- consensus_heatmap(res, k = 5)   

# hclust annotations
                                       
res_k5 <- res@object_list$`5`
cellpop <- rownames(res_k5$membership)
group_consensus <- rep(NA, length(cellpop))
names(group_consensus) <- cellpop
for (pop in cellpop) {
  group_consensus[pop] <- paste0("Eco", as.character(which(res_k5$membership[pop,]== max(res_k5$membership[pop,])))) }
names(group_consensus) <- c("Resident_macrophages_1", "CAF1", "EC_venous", "EC_arterial", "Resident_fibroblasts_1", "CAF3",                    
                          "TAM1", "Naive_Memory_T_cells", "TEC2", "Perivascular_TAM", "Exhausted_T_cells", "Inflammatory_macrophages",
                          "NK_like_T_cells", "GM1_ZG", "GM5_ECM", "GM2_ZF1", "GM4_ZR", "GM3_ZF2", "GM6_Translation", "GM7_Mitosis", "GM8_Hypoxia")
group_consensus2 <- group_consensus
group_consensus2[which(group_consensus == "Eco1")] <- "Eco2"
group_consensus2[which(group_consensus == "Eco2")] <- "Eco4"
group_consensus2[which(group_consensus == "Eco3")] <- "Eco5"
group_consensus2[which(group_consensus == "Eco4")] <- "Eco1"
group_consensus2[which(group_consensus == "Eco5")] <- "Eco3"
saveRDS(group_consensus2, "path/to/Ecotypes/results/Ecotypes_consensus.RDS")
                                        
cellpop <- p1$tree_col$labels
cellorder <- p1$tree_col$order
group_hclust <- cutree(p1$tree_col, k=5)
group_hclust2 <- group_hclust
group_hclust2[which(group_hclust == 1)] <- 4
group_hclust2[which(group_hclust == 2)] <- 5
group_hclust2[which(group_hclust == 4)] <- 1
group_hclust2[which(group_hclust == 5)] <- 2     

ngroups = 2
groups <- cutree(hclut_Test$tree_row, ngroups)
metadata <- cbind(dataclin_all[all_ACC, c("transcriptome", "cortisolsecretion", "OS_status", "cohort")], groups[all_ACC])
colnames(metadata) <- c("transcriptome_class", "cortisol_secretion", "death_event", "cohort", "Patient_cluster")
metadata$Patient_cluster[which(metadata$Patient_cluster == 1)] <- "C1A-like"
metadata$Patient_cluster[which(metadata$Patient_cluster == 2)] <- "C1B-like"
metadata$death_event <- as.character(metadata$death_event)
list_color <- list(transcriptome_class=c("C1A" = "darkmagenta", "C1B" ="hotpink3" ,"C2" = "deepskyblue4"),
                   cohort=c("TCGA" = "#386CB0", "FFPE" = "#BF5B17", "CIT" = "#7FC97F"),
                   death_event=c("0"="grey", "1"="black", "NA" = "white"),
                   cortisol_secretion=c("no"="grey", "yes"="black", "NA" = "white"),
                   Patient_cluster=c("C1A-like"="#440154FF", "C1B-like"="#238A8DFF"),
                  consensus_partition = c("Eco1"="#3f2d54", "Eco2"="#d18975", "Eco3"="#2d543d", "Eco4"="#8fd175", "Eco5"="#75b8d1"))                                       
                                        
# hclust plot
p_hclust <- pheatmap(Test, clustering_method = "ward.D2", 
                     cutree_rows = 2, cluster_cols = F,
                     annotation_row = metadata,  annotation_colors = list_color,
                     show_rownames = F, show_colnames = T, angle_col = 45)    

# stats                                      
fisher.test(metadata$Patient_cluster, metadata$transcriptome_class) # <10-16
fisher.test(metadata$Patient_cluster, metadata$cortisol_secretion) # <10-7
fisher.test(metadata$Patient_cluster, metadata$death_event) # <10-9
fisher.test(metadata$Patient_cluster, metadata$cohort)

# ecotypes scores
Eco2_ZF <- res_Z_all[, "GM2_ZF1"] + res_Z_all[, "GM3_ZF2"] + res_Z_all[, "GM1_ZG"] + res_Z_all[, "Exhausted_T_cells"] + res_Z_all[, "Resident_fibroblasts_1"] + res_Z_all[, "GM6_Translation"] 
Eco5_ZR <- res_Z_all[, "GM4_ZR"] + res_Z_all[, "Inflammatory_macrophages"] + res_Z_all[, "CAF1"]
Eco1_Hypoxia <- res_Z_all[, "GM8_Hypoxia"] + res_Z_all[, "GM7_Mitosis"] + res_Z_all[, "TEC2"] + res_Z_all[, "CAF3"]
Eco4_miscellanous2 <- res_Z_all[, "NK_like_T_cells"] + res_Z_all[, "Naive_Memory_T_cells"] + res_Z_all[, "Resident_macrophages_1"] + res_Z_all[, "TAM1"] + res_Z_all[, "EC_arterial"] + res_Z_all[, "GM5_ECM"]
Eco3_miscenallous1 <- res_Z_all[, "Perivascular_TAM"] + res_Z_all[, "EC_venous"]
res_Z_all <- cbind(res_Z_all, Eco1_Hypoxia, Eco2_ZF, Eco3_miscenallous1, Eco4_miscellanous2, Eco5_ZR) 


                                       
### Association of ecotypes with survival  -  Suppl Table 19 & 20

# univariate analysis
                                        
datatest <- cbind(annot[rownames(res_Z_all),], res_Z_all)                                          
Var <- c("stage", "grade", "cortisolsecretion", "Eco1_Hypoxia", "Eco2_ZF", "Eco3_miscenallous1" , "Eco4_miscellanous2", "Eco5_ZR")
ressurv <- matrix(nrow=0, ncol=9)
colnames(ressurv) <- c("variable", "HR_DFS", "CI_DFS", "p_DFS", "Cindex_DFS", "HR_OS", "CI_OS", "p_OS", "Cindex_OS")
                                        
for (v in Var) {
  coxDFS <- summary(coxph(Surv(datatest$follow_up_DFS_months, datatest$DFS_status) ~ datatest[ ,V]))
  coxOS <- summary(coxph(Surv(datatest$follow_up_OS_months, datatest$OS_status) ~ datatest[ ,V]))
  tmp <- matrix(nrow = nrow(coxOS$coefficients) + 1, ncol = 8)
  colnames(tmp) <- colnames(ressurv)
  tmp[1, "variable"] <- V
  tmp[2:nrow(tmp), "HR_DFS"] <- round(coxDFS$coefficients[,2],2)
  tmp[2:nrow(tmp), "CI_DFS"] <- paste0(round(coxDFS$conf.int[,3],2),"-",round(coxDFS$conf.int[,4],2))
  tmp[2:nrow(tmp), "p_DFS"] <- signif(coxDFS$coefficients[,5],3)
  tmp[2:nrow(tmp), "Cindex_DFS"] <- round(coxDFS$concordance[[1]], 3)
  tmp[2:nrow(tmp), "HR_OS"] <- round(coxOS$coefficients[,2],2)
  tmp[2:nrow(tmp), "CI_OS"] <- paste0(round(coxOS$conf.int[,3],2),"-",round(coxOS$conf.int[,4],2))
  tmp[2:nrow(tmp), "p_OS"] <- signif(coxOS$coefficients[,5],3)
  tmp[2:nrow(tmp), "Cindex_OS"] <- round(coxOS$concordance[[1]], 3)
  ressurv <- rbind(ressurv, tmp) }                                        
                                        
# multivariate analysis
                                        
complete_model <- coxph(Surv(follow_up_OS_months, OS_status) ~ Eco1_Hypoxia + Eco2_ZF + Eco5_ZR + stage, data = datatest)
step_model_OS <- stepAIC(complete_model_OS, direction = "both")
summary(step_model_OS)                                        

complete_model_DFS <- coxph(Surv(follow_up_DFS_months, DFS_status) ~ Eco1_Hypoxia + Eco2_ZF + Eco5_ZR + stage, data = datatest)
step_model_DFS <- stepAIC(complete_model_DFS, direction = "both")
summary(step_model_DFS)    

                                        
### Best survival model - comparison of C indexes & nested models - Suppl Table 20

resCindex <- matrix(nrow=0, ncol=3)                                        
colnames(resCindex) <- c("variable", "Cindex_DFS", "Cindex_OS")
                                        
# stage + grade

coxDFS <- summary(coxph(Surv(follow_up_DFS_months, DFS_status) ~ stage + grade, data = datatest))
coxOS <- summary(coxph(Surv(follow_up_OS_months, OS_status) ~ stage + grade, data = datatest))
tmp <- c("stage + grade", round(coxDFS$concordance[[1]], 3), round(coxOS$concordance[[1]], 3))
resCindex <- rbind(resCindex, tmp)

#  stage + grade + secretion

coxDFS <- summary(coxph(Surv(follow_up_DFS_months, DFS_status) ~ stage + grade + cortisolsecretion, data = datatest))
coxOS <- summary(coxph(Surv(follow_up_OS_months, OS_status) ~ stage + grade + cortisolsecretion, data = datatest))
tmp <- c("stage + grade + cortisol", round(coxDFS$concordance[[1]], 3), round(coxOS$concordance[[1]], 3))
resCindex <- rbind(resCindex, tmp)

# ecotypes
coxDFS <- summary(coxph(Surv(follow_up_DFS_months, DFS_status) ~ Eco1_Hypoxia + Eco2_ZF + Eco5_ZR, data = datatest))
coxOS <- summary(coxph(Surv(follow_up_OS_months, OS_status) ~ Eco1_Hypoxia + Eco2_ZF + Eco5_ZR, data = datatest))
tmp <- c("Eco1 + Eco2 + Eco5", round(coxDFS$concordance[[1]], 3), round(coxOS$concordance[[1]], 3))
resCindex <- rbind(resCindex, tmp)

# combined model
coxDFS <- summary(coxph(Surv(follow_up_DFS_months, DFS_status) ~ stage + grade + cortisolsecretion + Eco1_Hypoxia + Eco2_ZF + Eco5_ZR, data = datatest))
coxOS <- summary(coxph(Surv(follow_up_OS_months, OS_status) ~ stage + grade + cortisolsecretion + Eco1_Hypoxia + Eco2_ZF + Eco5_ZR, data = datatest))
tmp <- c("stage + cortisol + Eco1 + Eco2 + Eco5", round(coxDFS$concordance[[1]], 3), round(coxOS$concordance[[1]], 3))
resCindex <- rbind(resCindex, tmp)
                                                                  
# nested models DFS
cox_base <- coxph(Surv(follow_up_DFS_months, DFS_status) ~ stage + grade + cortisolsecretion, data = datatest)
cox_complete <- coxph(Surv(follow_up_DFS_months, DFS_status) ~ stage + grade + cortisolsecretion + Eco1_Hypoxia + Eco2_ZF + Eco5_ZR, data = datatest)
lrtest(cox_complete, cox_base)

# nested models OS
cox_base <- coxph(Surv(follow_up_OS_months, OS_status) ~ stage + grade + cortisolsecretion, data = datatest)
cox_complete <- coxph(Surv(follow_up_OS_months, OS_status) ~ stage + grade + cortisolsecretion + Eco1_Hypoxia + Eco2_ZF + Eco5_ZR, data = datatest)
lrtest(cox_complete, cox_base)                                      
                                        
