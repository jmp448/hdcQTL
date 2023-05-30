import scanpy as sc
import pandas as pd
import numpy as np
import re

lineage_adata = snakemake.input['adata']
pseudotime_adata = snakemake.output['adata_pseudotime']
nbins = int(re.findall(r'\d+', str(snakemake.wildcards['nbins']))[0])

adata = sc.read_h5ad(lineage_adata)
sc.tl.pca(adata)

bin_names = ['bin' + str(i) for i in range(1, nbins + 1)]
adata.obs['pseudotime'] = adata.obsm['X_pca'][:, 0]
adata.obs['stage'] = pd.cut(adata.obs['pseudotime'], nbins, labels=bin_names)

# Trim the bins at either end, reassigning all cells to the nearest bin
adata.obs.loc[adata.obs['stage'] == bin_names[0], 'stage'] = bin_names[1]
adata.obs.loc[adata.obs['stage'] == bin_names[-1], 'stage'] = bin_names[-2]

adata.write_h5ad(pseudotime_adata)
