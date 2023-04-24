import scanpy as sc
import pandas as pd
import numpy as np

raw_counts_loc = snakemake.input[0]
normalized_adata_all_loc = snakemake.output[0]
normalized_adata_hvg_loc = snakemake.output[1]

# Load data
adata = sc.read_h5ad(raw_counts_loc)
adata.layers['counts'] = adata.X.copy()

# Log normalize
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers['log1pPF'] = adata.X.copy()

# Apply additional proportional fitting
sc.pp.normalize_total(adata)
adata.layers['PFlog1pPF'] = adata.X.copy()

# Identify highly variable genes
sc.pp.highly_variable_genes(adata, layer='log1pPF', n_top_genes=5000)

# Save full object
adata.write_h5ad(normalized_adata_all_loc)

# Filter to highly variable genes and save a lighter-weight object that's subset to them
adata_hvg = adata[:, adata.var['highly_variable']]
adata_hvg.write_h5ad(normalized_adata_hvg_loc)
