import scanpy as sc
import pandas as pd
import numpy as np

lineage_adata = snakemake.input['adata']
stage_labels = snakemake.input['stage_labels']
pseudotime_adata = snakemake.output['adata_pseudotime']

adata = sc.read_h5ad(lineage_adata)
sc.tl.pca(adata)

stage_labels = pd.read_csv(stage_labels)['stages']
adata.obs['pseudotime'] = adata.obsm['X_pca'][:, 0]
adata.obs['stage'] = pd.cut(adata.obs['pseudotime'], 5, labels=stage_labels)

adata.write_h5ad(pseudotime_adata)
