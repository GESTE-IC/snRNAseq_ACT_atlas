library(Seurat)
library(ggplot2)
library(patchwork)
library(viridis)

ACC2a <- readRDS("/path/to/results/ACC2a_annot.RDS")

figure_dir <- "/path/to/figures/"
sample_name <- "ACC2a_"

SpatialDimPlot(ACC2a, alpha = 0, crop = FALSE) & 
  ggtitle("") &
  NoLegend()

SpatialDimPlot(ACC2a, group.by = "SCT_snn_res.0.3", crop = FALSE, pt.size.factor = 1.3, image.alpha = 0) &
  scale_fill_viridis_d(option = "viridis", direction = -1) & 
  ggtitle("") &
  theme(legend.key.size = unit(2, "cm"),
        legend.title = element_text(size = 30),
        legend.text = element_text(size = 30), 
        legend.key = element_blank()) &
  guides(fill = guide_legend(override.aes = list(size = 20))) &
  labs(fill = "Clusters res 0.3")

SpatialFeaturePlot(ACC2a, features = "SPP1", image.alpha = 0, crop = FALSE, pt.size.factor = 1.3) &
  scale_fill_viridis_c(option = "plasma")

SpatialFeaturePlot(ACC2a, features = "CYP11A1", image.alpha = 0, crop = FALSE, pt.size.factor = 1.3) &
  scale_fill_viridis_c(option = "plasma")

SpatialFeaturePlot(ACC2a, features = c("CAF1", "TAM1", "NK.like.T.cells", "Inflammatory.macrophages"), ncol = 2, crop = FALSE, pt.size.factor = 1.3, image.alpha = 0) & 
  scale_fill_viridis_c(option = "turbo") &
  theme(legend.key.size = unit(3, "cm"),
        legend.title = element_text(size=30),
        legend.text = element_text(size=30), 
        legend.key = element_blank())

SpatialFeaturePlot(ACC2a, features = c("GM1_ZG_1", "GM2_ZF1_2", "GM3_ZF2_3", "GM4_ZR_4", "GM5_ECM_5", "GM6_Translation_6", "GM7_Mitosis_7", "GM8_Hypoxia_8"), ncol = 4, crop = FALSE, pt.size.factor = 1.3, image.alpha = 0) &
  scale_fill_viridis_c(option = "turbo") &
  theme(legend.key.size = unit(3, "cm"),
        legend.title = element_text(size=30),
        legend.text = element_text(size=30), 
        legend.key = element_blank())
