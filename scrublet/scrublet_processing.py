# Scrublet v0.2.3

import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import os
import numba
import pandas as pd
from os import walk

plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rc('font', size=14)
plt.rcParams['pdf.fonttype'] = 42

input_path = "/path/to/files.mtx"

from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(input_path) if isfile(join(input_path, f))]

for Raw_matrix in onlyfiles:
    output_path = "/path/to/results/Scrublet/figures/"
    print(Raw_matrix)
    counts_matrix = scipy.io.mmread(input_dir + '/' + Raw_matrix).T.tocsc()
    print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
    
    scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
    scrub.plot_histogram();
    plt.savefig(output_path+'Hist_'+Raw_matrix.strip('.mtx')+'.png')

    print('Running UMAP...')
    scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
    print('Done.')
    scrub.plot_embedding('UMAP', order_points=True);4
    plt.savefig(output_path+'Umap_'+Raw_matrix.strip('.mtx')+'.png')
    
    print('{} doublets'.format(np.sum(predicted_doublets)))
    print('{} doublets indexes'.format([i for i, x in enumerate(predicted_doublets) if x]))

    df = pd.DataFrame({
        'doublet_score': scrub.doublet_scores_obs_,
        'predicted_doublet': scrub.predicted_doublets_
    })
    output_dir = "/path/to/results/Scrublet/outputs/"
    df.to_csv(output_path+Raw_matrix.strip('.mtx')+'_table.csv', index=False)
    
