import scanpy as sc
import numpy as np
import pandas as pd
import anndata
from sklearn.preprocessing import OneHotEncoder
import random

random.seed(1234)

'''
This script is used to generate pseudobulk cells from single cells EB expression data. 
The input are (1) anndata object with normalized counts and (2) a matched scVI embedding file. 
The output are: 
(1) cluster_asignments (dataframe of mapping between cell id and cluster assignment) that are 
used in generate pseudobulk expression.
(2) anndata object for the pseudobulk expressions
'''

adata_raw_loc = snakemake.input['adata_raw']
pseudocell_mapping_loc = snakemake.input['pseudocell_mapping']
adata_pseudocell_loc = snakemake.output['adata_pseudocell']
pseudocell_names = snakemake.output['adata_pseudocell_names']
gene_names = snakemake.output['adata_gene_names']

adata = sc.read_h5ad(adata_raw_loc)
pseudocell_mapping = pd.read_csv(pseudocell_mapping_loc, sep="\t")

adata = adata[pseudocell_mapping['cell'], :] 
adata.obs['pseudocell'] = list(pseudocell_mapping['pseudocell_15'])
adata.obs['pseudocell'] = adata.obs['pseudocell'].astype("category")

# Aggregate (sum) across cells within a pseudocell
onehot = OneHotEncoder(sparse=True).fit_transform(adata.obs[['pseudocell']])
pseudocell_exp = onehot.transpose() * adata.X # (adata.X.transpose() * onehot).transpose

# Generate obs and var for new anndata object
## gene features
pseudocell_var = adata.var

## sample (pseudo cell) features
pseudocell_obs = pd.DataFrame(pseudocell_mapping[['pseudocell_15']].value_counts()).reset_index().rename(columns={'pseudocell_15': 'pseudocell', 0: 'ncells'}) 
pseudocell_obs['donor_id'] = [s.split("_")[0] for s in pseudocell_obs['pseudocell']]
pseudocell_obs = pseudocell_obs.set_index("pseudocell")
pseudocell_obs = pseudocell_obs.reindex(adata.obs['pseudocell'].astype("category").cat.categories)

# make sure indices match
assert np.array_equal(pseudocell_obs.index, adata.obs['pseudocell'].astype("category").cat.categories)

# make an anndata object
adata_pseudocell = anndata.AnnData(X=pseudocell_exp, dtype=int, obs=pseudocell_obs, var=pseudocell_var)
sc.pp.filter_genes(adata_pseudocell, min_cells=10)

adata_pseudocell.write_h5ad(adata_pseudocell_loc)

# need to save cell and gene names to a tsv since they get disrupted during file conversion
pd.DataFrame(adata_pseudocell.obs_names).to_csv(pseudocell_names, sep='\t', header=False, index=False)
pd.DataFrame(adata_pseudocell.var_names).to_csv(gene_names, sep='\t', header=False, index=False)
