import os

import numpy as np
import scanpy as sc

# check if file is downloaded
# https://ncbi.nlm.nih.gov/geo/download/?acc=GSM3852755&format=file&file=GSM3852755%5FE15%5F5%5Fcounts%2Etar%2Egz
if not os.path.exists("GSM3852755_E15_5_counts.tar.gz"):
    print("File not found. Please download the file.")
    quit()

# set seed to 5005
np.random.seed(5005)

adata = sc.read_10x_mtx("GSM3852755_E15_5_counts.tar.gz")

# quality conrol
# filter low quality genes and cells
# cells that do not have many genes expressed
# geens that are not expressed in many cells

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# calculate metrics to detemine low quality and doublets
# total genes counts
# percent of gene counts from mitochondrial genes

adata.obs["total_counts"] = adata.X.sum(axis=1)
adata.obs["pct_counts_mt"] = np.sum(adata[:, adata.var_names.str.startswith("MT-")].X, axis=1) / adata.obs["total_counts"]

# normalize and log+1 tranformation of data
# determine highly variable genesy

adata = sc.pp.normalize_total(adata, target_sum=1e4)
adata = sc.pp.log1p(adata)
adata = sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)


# running PCA, UMAP, and TSNE for viualization
# determine the appropraite number of princial components with an elbow plot

adata = sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")
sc.pl.pca_variance_ratio(adata, log=True)


# calculate neighborhood for 2D visualization
# leiden clustering to find cell types
# use resolution parameters to detemine number of clusters

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color=["leiden", "pct_counts_mt"], legend_loc="on data")

# calculate which genes define each cluster
# do these genes efine cell types in literature? eg GSEA, initial paper

sc.tl.score_genes(adata, "leiden", "leiden_genes", gene_list=adata.var_names)