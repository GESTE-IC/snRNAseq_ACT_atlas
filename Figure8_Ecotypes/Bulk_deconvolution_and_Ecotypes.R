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

res_cibersort_TCGA <-  read.csv2("/path/to/Cibersort/results/TCGA/CIBERSORTx_Adjusted.txt", sep="\t")
rownames(res_cibersort_TCGA) <- res_cibersort_TCGA$Mixture
celltypes <- setdiff(colnames(res_cibersort_TCGA), c("Mixture","P.value","Correlation","RMSE"))
res_cibersort_TCGA[,celltypes] <- apply(res_cibersort_TCGA[,celltypes], 2, function(x) as.numeric(x))
                                       
res_cibersort_CIT <-  read.csv2("/path/to/Cibersort/results/CIT/CIBERSORTx_Adjusted.txt", sep="\t")
rownames(res_cibersort_CIT) <- res_cibersort_CIT$Mixture
res_cibersort_CIT[,celltypes] <- apply(res_cibersort_CIT[,celltypes], 2, function(x) as.numeric(x))
                                        
res_cibersort_FFPE <-  read.csv2("/path/to/Cibersort/results/FFPE/CIBERSORTx_Adjusted.txt", sep="\t")
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

                              
### Association of single cell signatures with survival  -  Suppl Table 7  

# univariate analyses                                        
Varqual <- c("age", "sex", "cortisolsecretion", "stade", "mitoses", "ki67")
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
  tmp[1, "N"] <- NpourV
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
                          GM3_ZF2 + GM7_Mitosis + stade, 
                        data = datatest)
step_model_OS <- stepAIC(complete_model_OS, direction = "both")                                        

complete_model_DFS <- coxph(Surv(follow_up_DFS_months, DFS_status) ~ CAF1 + EC_arterial + Resident_fibroblasts_1 +
                          CAF3 + TEC2 + Inflammatory_macrophages + 
                          GM3_ZF2 + GM7_Mitosis + stade, 
                        data = datatest)
step_model_DFS <- stepAIC(complete_model_DFS, direction = "both")                                 
                                        
