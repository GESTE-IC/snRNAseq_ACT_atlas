library(GSVA)

list_GM <- readRDS("/path/to/GeneModules/results/Gene_modules.RDS")

### load dataset bulk
counts_FFPE <- read.table("/path/to/public_data/counts_FFPE.txt"), sep="\t", header=T)
annot_FFPE <- read.csv2("/path/to/public_data/annot_FFPE.csv"), dec=".")
counts_TCGA <- read.table("/path/to/public_data/counts_TCGA.txt"), sep="\t", header=T)
annot_TCGA <- read.csv2("/path/to/public_data/annot_TCGA.csv"), dec=".")
counts_CIT <- read.table("/path/to/public_data/counts_CIT.txt"), sep="\t", header=T)
annot_CIT <- read.csv2("/path/to/public_data/annot_CIT.csv"), dec=".")

### ssgsea scores

res_deconv <- list()
for (i in 1:8) {
  param_FFPE <- ssgseaParam(counts_FFPE, list_prog[i])
  res_deconv[["FFPE"]] <- gsva(param_FFPE)
  param_TCGA <- ssgseaParam(counts_TCGA, list_prog[i])
  res_deconv[["TCGA"]] <- gsva(param_TCGA)
  param_CIT <- ssgseaParam(counts_CIT, list_prog[i])
  res_deconv[["CIT"]] <- gsva(param_CIT) }

saveRDS(res_deconv, ""/path/to/Deconvolution/results/Deconvolution_GM.RDS")
