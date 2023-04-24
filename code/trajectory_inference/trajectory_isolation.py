import scanpy as sc
import numpy as np
import pandas as pd
import anndata
from statsmodels.stats.multitest import multipletests
import scdrs

adata_loc = snakemake.input['adata']
gs_loc = snakemake.input['gs']
scdrs_output_loc = snakemake.output['scdrs']
filtered_anndata_loc = snakemake.output['adata_filtered']

# Load and preprocess EB anndata
eb_adata = sc.read_h5ad(adata_loc)

# Load gene set
gs = pd.read_csv(gs_loc, sep="\t")
markers = gs['Gene']
weights = gs['Weight']

# Score cells with scDRS
scores = scdrs.score_cell(eb_adata, gene_list=markers, gene_weight=weights, n_ctrl=20)
scores['pval_adj_bh'] = multipletests(scores["pval"].values, method="fdr_bh")[1]
scores.to_csv(scdrs_output_loc)

# Filter the anndata to cells along this trajectory
eb_adata.obs['trajectory_inclusion_fdr'] = scores['pval_adj_bh']
eb_adata = eb_adata[eb_adata.obs['trajectory_inclusion_fdr'] <= 0.1] 
eb_adata.write_h5ad(filtered_anndata_loc)
