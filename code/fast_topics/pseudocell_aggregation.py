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
adata.obs['pseudocell'] = list(pseudocell_mapping['individual_leiden_clusters_15'])
adata.obs['pseudocell'] = adata.obs['pseudocell'].astype("category")

# Aggregate (sum) across cells within a pseudocell
onehot = OneHotEncoder(sparse=True).fit_transform(adata.obs[['pseudocell']])
pseudocell_exp = onehot.transpose() * adata.X # (adata.X.transpose() * onehot).transpose

# Generate obs and var for new anndata object
## gene features
pseudocell_var = adata.var

## sample (pseudo cell) features
pseudocell_obs = pd.DataFrame(pseudocell_mapping[['individual_leiden_clusters_15']].value_counts()).reset_index().rename(columns={'individual_leiden_clusters_15': 'pseudocell', 0: 'ncells'}) 
pseudocell_obs['donor_id'] = [s.split("_")[0] for s in pseudocell_obs['pseudocell']]
pseudocell_obs = pseudocell_obs.set_index("pseudocell")
pseudocell_obs = pseudocell_obs.reindex(adata.obs['pseudocell'].astype("category").cat.categories)

# make sure indices match
assert np.array_equal(pseudocell_obs.index, adata.obs['pseudocell'].astype("category").cat.categories)

# make an anndata object
adata_pseudocell = anndata.AnnData(X=pseudocell_exp, dtype=int, obs=pseudocell_obs, var=pseudocell_var)
                      
adata_pseudocell.write_h5ad(adata_pseudocell_loc)

# need to save cell and gene names to a tsv since they get disrupted during file conversion
pd.DataFrame(adata_pseudocell.obs_names).to_csv(pseudocell_names, sep='\t', header=False, index=False)
pd.DataFrame(adata_pseudocell.var_names).to_csv(gene_names, sep='\t', header=False, index=False)

# def generate_sum_cluster_pseudobulk_expression(adata, ordered_pseudobulk_samples, cluster_assignments):
#     # the adata should be counts data
#     num_samples = len(ordered_pseudobulk_samples)
#     num_genes = adata.shape[1]
#     pseudobulk_expr = np.zeros((num_samples, num_genes))
# 
#     for pseudobulk_sample_num, pseudobulk_sample_name in enumerate(ordered_pseudobulk_samples):
#         # Get cell indices corresponding to this pseudobulk sample
#         indices = cluster_assignments == pseudobulk_sample_name
#         # Fill in pseudobulk expr
#         pseudobulk_expr[pseudobulk_sample_num, :] = np.asarray(np.sum(adata.X[indices,:], axis=0))
#     return pseudobulk_expr
# 
# 
# # dataframe of mapping between cell id and cluster assignment
# cluster_asignments = adata_woPilot_sampled.obs[["individual_leiden_clusters_15"]]
# 
# # order pseudocell ids
# ordered_pseudobulk_samples_raw = sorted(np.unique(adata_woPilot_sampled.obs.individual_leiden_clusters_15),
#                            key= lambda i: (int(i.lstrip('NA').split(':')[0]), int(i.split(':')[1])))
# 
# # save the cell-pseudocell mapping
# cluster_asignments.to_csv(snakemake_output[0],sep='\t')
# 
# 
# def log1pnormalized2counts(adata, scale_by='sizeFactor'):
#     adata.X.data = np.exp(adata.X.data) - 1  # undo log1
#     inplace_row_scale(adata.X, adata.obs[scale_by].values)  # tocounts
#     adata.X.data = adata.X.data.astype(int)
# 
# '''
# Generate pseudobulk level expression anndata
# '''
# # Convert normalized counts to raw counts
# log1pnormalized2counts(adata_woPilot_sampled, scale_by='sizeFactor')
# 
# # # make this dataframe a series for generate_sum_cluster_pseudobulk_expression function parameter purpose
# cluster_asignments = cluster_asignments["individual_leiden_clusters_15"]
# pseudo_expr = generate_sum_cluster_pseudobulk_expression(adata_woPilot_sampled, ordered_pseudobulk_samples_raw, cluster_asignments)
# 
# # generate the gene feature and sample information
# ## gene features
# var = adata_woPilot_sampled.var
# 
# ## sample (pseudo cell) features
# donor_id = [i.split(":")[0] for i in ordered_pseudobulk_samples_raw]
# sample_feature_dic = {"pseudo_cell":ordered_pseudobulk_samples_raw,
#                      "donor_id": donor_id}
# obs = pd.DataFrame(sample_feature_dic)
# obs = obs.set_index("pseudo_cell")
# 
# # make an anndata object
# adata_woPilot_sampled_pseudobulk = anndata.AnnData(X=pseudo_expr,dtype=int, 
#                       obs=obs, 
#                       var=var)
# 
# # add ncells per pseudocell in the obs feature
# adata_woPilot_sampled_pseudobulk.obs["ncells"] = pd.Series(cluster_asignments.value_counts())
# 
# # save the anndata object
# adata_woPilot_sampled_pseudobulk.write(snakemake_output[1])
