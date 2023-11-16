#!/usr/bin/env python3.9

import sys
import numpy as np
import pandas as pd
import scanpy as sc
import os
import anndata as ad

## Output versions
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

### Load Data of 5 reference donors
src = sys.argv[1]
dest = sys.argv[2]

for filename in os.listdir(src):

    data = sc.read_h5ad(src+"/"+filename)

    ## Preprocessing data

    # Makes the index unique by appending a number string to each duplicate index element: ‘1’, ‘2’, etc.
    data.var_names_make_unique()
    # Plot high expression genes
    sc.pl.highest_expr_genes(data, n_top=50,)
    # Filter low expression genes and cells
    sc.pp.filter_cells(data, min_genes=50)
    #sc.pp.filter_cells(data, min_cells=3)
    # Check doublets
    sc.external.pp.scrublet(data)
    # Delete doublets
    data = data[data.obs['predicted_doublet'] == False] 

    # Filter based on mitochondrial counts and total counts
    data.var['mt'] = data.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    # Visualize qc
    import matplotlib.pyplot as plt
    mito_filter = 3
    upper_lim = np.quantile(data.obs.n_genes_by_counts.values, .98)
    lower_lim = np.quantile(data.obs.n_genes_by_counts.values, .02)
    fig, axs = plt.subplots(ncols = 2, figsize = (8,4))
    #sc.pl.scatter(data, x='total_counts', y='pct_counts_mt', ax = axs[0], show=False)
    sc.pl.scatter(data, x='total_counts', y='n_genes_by_counts', ax = axs[1], show = False)
    #draw horizontal red lines indicating thresholds.
    axs[0].hlines(y = mito_filter, xmin = 0, xmax = max(data.obs['total_counts']), color = 'red', ls = 'dashed')
    axs[1].hlines(y = upper_lim, xmin = 0, xmax = max(data.obs['total_counts']), color = 'red', ls = 'dashed')
    axs[1].hlines(y = lower_lim, xmin = 0, xmax = max(data.obs['total_counts']), color = 'red', ls = 'dashed')
    fig.tight_layout()
    plt.show()

    # Delete outsider counts
    data = data[(data.obs.n_genes_by_counts < upper_lim) & (data.obs.n_genes_by_counts > lower_lim)]
    #data = data[data.obs.pct_counts_mt < mito_filter, :]


    ## Normalization
    # Convert Scipy to array
    data.X=data.X.A
    # Normalize every cell to 10,000 UMI and then change to log counts
    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    # Check data
    print(data.X.sum(axis=1))
    # Store precessed data
    data.raw = data

    ### Clustering
    # Regress out unwanted sources of variation
    sc.pp.regress_out(data, ['total_counts', 'pct_counts_mt'])
    print(data)
    # Normalize each gene to the unit variance of that gene
    sc.pp.scale(data,max_value=10)
    # Calculate PCA (default is 50 pcs)
    sc.tl.pca(data, svd_solver='arpack')
    # Elbow plot to show the pca significance
    sc.pl.pca_variance_ratio(data, log=True, n_pcs=50)
    plt.savefig(dest+'/Elbow_Plot_for_PCA_Significance.png')
    plt.close()
    # From around 30th PC, there is no significant difference, so take pc=30 for neighbors computing
    sc.pp.neighbors(data, n_pcs = 30)
    # Visulize the data with umap
    sc.tl.umap(data)
    sc.pl.umap(data)
    # Clustering with Leiden algorithm (resolution should be adjusted later)
    sc.tl.leiden(data, resolution = 0.5)
    sc.pl.umap(data, color=['leiden'])

    ### Annotation
    # Get the marker genes based on Leiden
    sc.tl.rank_genes_groups(data, 'leiden')
    # Plot the marker genes
    sc.pl.rank_genes_groups(data, n_genes=20, sharey=False)
    plt.savefig(dest+'/Rank_Genes.png')
    plt.close()
    # Save the Marker genes
    markers_data = sc.get.rank_genes_groups_df(data, None)
    # Filter marker genes with p-value
    markers_data = markers_data[(markers_data.pvals_adj < 0.05)&(markers_data.logfoldchanges > .5)]
    print(markers_data)

    # Plot umap with legend on it   
    sc.pl.umap(data, color = ['leiden'], frameon = False, legend_loc = "on data", title="Leiden Clustering umap")
    plt.savefig(dest+'/Leiden_umap.png')
    plt.close()



