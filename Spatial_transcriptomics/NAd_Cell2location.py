#!/usr/bin/env python3

import sys 
import scanpy as sc 
import anndata 
import pandas as pd 
import numpy as np 
import os 
import gc

os.environ["THEANO_FLAGS"] = 'device=cuda0,floatX=float32,force_device=True'
# if using the CPU uncomment this:
#os.environ["THEANO_FLAGS"] = 'device=cpu,floatX=float32,openmp=True,force_device=True'

import cell2location
import squidpy as sq

import matplotlib as mpl
from matplotlib import rcParams
import matplotlib.pyplot as plt
import seaborn as sns

st_data_folder = "/path/to/NAD4/count/matrix"
sc_data_folder = "/path/to/cell2location_data"
results_data_folder = "/path/to/cell2location_results"
ref_run_name = f'{results_data_folder}/reference_signatures'
run_name = f'{results_data_folder}/cell2location_map'

adata_vis = sq.read.visium(path=st_data_folder, counts_file="NAD4_filtered_feature_bc_matrix.h5", source_image_path=st_data_folder + "/spatial/tissue_hires_image.png")
adata_vis.var_names_make_unique()
adata_vis.var['SYMBOL'] = adata_vis.var_names

# import ref data
adata_ref = sc.read_csv(sc_data_folder + "/singlecell_ref_count_matrix.tsv", delimiter="\t", first_column_names=True)
ref_metadata = pd.read_csv(sc_data_folder + "/singlecell_ref_annotation_table.tsv", sep="\t")
adata_ref.obs["bio_celltype"] = pd.Categorical(ref_metadata["bio_celltype"])

# filtering of important genes
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

# data preparation : identification of covariates for regression model 
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        # batch_key='Sample',
                        # cell type, covariate used for constructing signatures
                        labels_key="bio_celltype",
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        # categorical_covariate_keys=['donor_name']
                       )

from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
mod.view_anndata_setup()

# You are using a CUDA device ('NVIDIA GeForce RTX 3070') that has Tensor Cores. To properly utilize them, you should set `torch.set_float32_matmul_precision('medium' | 'high')` which will trade-off precision for performance. For more details, read https://pytorch.org/docs/stable/generated/torch.set_float32_matmul_precision.html#torch.set_float32_matmul_precision
from torch import set_float32_matmul_precision
set_float32_matmul_precision('high')
mod.train(max_epochs=250, use_gpu=True)
mod.plot_history(20)
plt.show()
plt.savefig("loess_curve_reg_model.png")

# Export the summary of the posterior distrib
adata_ref = mod.export_posterior(adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True})

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)
mod.plot_QC()

# extracxt reference cell types as a pd.DataFrame
# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, :]

# SPATIAL MAPPING
# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis)

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=5,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=20
)
mod.view_anndata_setup()

mod.train(max_epochs=30000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True
         )
# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);
plt.show()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

# Save model
mod.save(f"{run_name}", overwrite=True)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file

# Load model when needed
adata_file = f"{run_name}/sp.h5ad"
adata_vis = sc.read_h5ad(adata_file)
mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# save matrix
abundance_file = f"{run_name}/abundance_matrix.csv"
ab_mat = pd.DataFrame(adata_vis.obsm["q05_cell_abundance_w_sf"])
ab_mat["sample_name"] = ab_mat.index
ab_mat.to_csv(abundance_file, index=False)
