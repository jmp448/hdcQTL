import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import scvi
import math
import matplotlib.pyplot as plt
import scipy.sparse
from glob import glob
import re
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LinearRegression

# raw_loc = "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/highpass/eb_raw.qc.h5ad"
# pseudotime_annotated = "/project2/gilad/jpopp/ebQTL/data/trajectory_inference/cm_lineage/eb_cm_lineage.pseudotime.adata"
# 
# celltype_summary_loc = "/project2/gilad/jpopp/ebQTL/data/static_qtl_calling/eb_cmstages/pseudobulk_tmm/samples_per_celltype.tsv"
# sample_summary_loc = "/project2/gilad/jpopp/ebQTL/data/static_qtl_calling/eb_cmstages/pseudobulk_tmm/sample_summary.tsv"
# pseudobulk_loc = "/project2/gilad/jpopp/ebQTL/data/static_qtl_calling/eb_cmstages/pseudobulk_tmm/eb_cmstages.pseudobulk_tmm.tsv"

raw_loc = snakemake.input['raw']
pseudotime_annotated = snakemake.input['pseudotimed']

donor_summary_loc = snakemake.output['donor_summary']
sample_summary_loc = snakemake.output['sample_summary']
pseudobulk_loc = snakemake.output['pseudobulk']

nbins = int(re.findall(r'\d+', str(snakemake.wildcards['nbins']))[0])
is_trimmed = str(snakemake.wildcards['nbins'])[-7:]=='trimmed'
drops_allowed = int(snakemake.params['drops_allowed'])
min_cells = int(snakemake.params['min_cells_per_bin'])

# Load data (then dump it to avoid overusing memory)
adata_full = sc.read_h5ad(raw_loc)
adata_lineage = sc.read_h5ad(pseudotime_annotated)

adata = adata_full[adata_lineage.obs_names]
adata.obs['stage'] = adata_lineage.obs['stage']
adata.obs['pseudotime'] = adata_lineage.obs['pseudotime']

del adata_full
del adata_lineage

## Create summary table
cell_counts = adata.obs[['donor_id', 'stage']]
cell_counts = pd.DataFrame(cell_counts.groupby('stage').value_counts()).reset_index(inplace=False).rename(columns={0: "n_cells_unfiltered", 'stage': 'type'})

cell_counts['individual'] = [s.replace("NA", "") for s in cell_counts['donor_id']]
cell_counts['ind_type'] = cell_counts['individual'].astype(str) + "_" + cell_counts['type'].astype(str)
cell_counts = cell_counts[['ind_type', 'individual', 'type', 'n_cells_unfiltered']]

# drop samples w less than 5 cells
cell_counts['dropped'] = cell_counts['n_cells_unfiltered'] < min_cells
samples_inc_1 = cell_counts[(cell_counts['dropped'] == False)]['ind_type']

cell_subset_1 = adata.obs[['donor_id', 'stage']].copy().rename(columns={'stage': 'type'})
cell_subset_1['ind'] = [s.replace("NA", "") for s in cell_subset_1['donor_id'].astype(str)]
cell_subset_1['sample'] = cell_subset_1['ind'].astype(str) + "_" + cell_subset_1['type'].astype(str)
cell_subset_1 = cell_subset_1[cell_subset_1['sample'].isin(samples_inc_1)]

adata = adata[cell_subset_1.index]

## Update summary table
filtered_counts = adata.obs[['donor_id', 'stage', 'total_counts', 'pseudotime']].copy()
filtered_counts['n_cells_filtered'] = 1
filtered_counts['individual'] = [s.replace("NA", "") for s in filtered_counts['donor_id']]
filtered_counts['ind_type'] = filtered_counts['individual'].astype(str) + "_" + filtered_counts['stage'].astype(str)
filtered_counts = filtered_counts.drop(columns=['donor_id', 'individual', 'stage'])
filtered_counts = filtered_counts.groupby('ind_type').agg({'total_counts': 'sum', 'n_cells_filtered': 'count', 'pseudotime': 'median'})
filtered_counts = filtered_counts.reset_index().astype({'total_counts': 'int'})

cell_counts_filtered = cell_counts.merge(filtered_counts, on='ind_type', how='left').fillna({'total_counts': 0, 'n_cells_filtered': 0}).astype({'total_counts': 'int', 'n_cells_filtered': 'int'})
cell_counts_filtered['dropped'] = cell_counts_filtered['n_cells_filtered'] < min_cells
cell_counts_filtered = cell_counts_filtered.sort_values(by="n_cells_filtered", ascending=False)

# drop donors that are missing in more than the cutoff number of pseudotime bins
bin_counts_filtered = cell_counts_filtered[cell_counts_filtered['n_cells_filtered'] >= min_cells]
bin_counts_filtered = pd.DataFrame(bin_counts_filtered[['individual']].value_counts()).reset_index(inplace=False).rename(columns={0: "n_filtered"})
if is_trimmed:
  bin_counts_filtered = bin_counts_filtered[bin_counts_filtered['n_filtered'] >= nbins-drops_allowed-2]  
else:
  bin_counts_filtered = bin_counts_filtered[bin_counts_filtered['n_filtered'] >= nbins-drops_allowed] 
bin_counts_filtered.to_csv(donor_summary_loc, sep="\t", index=False)

donors_inc = bin_counts_filtered['individual']
cell_counts_filtered['dropped'] = (cell_counts_filtered['dropped']) | ~(cell_counts_filtered['individual'].isin(donors_inc))
samples_inc_2 = cell_counts_filtered[(cell_counts_filtered['dropped'] == False)]['ind_type']
cell_counts_filtered.to_csv(sample_summary_loc, sep="\t")

## Aggregation
cell_subset_2 = adata.obs[['donor_id', 'stage']].copy().rename(columns={'stage': 'type'})
cell_subset_2['ind'] = [s.replace("NA", "") for s in cell_subset_2['donor_id'].astype(str)]
cell_subset_2['sample'] = cell_subset_2['ind'].astype(str) + "_" + cell_subset_2['type'].astype(str)
cell_subset_2 = cell_subset_2[cell_subset_2['sample'].isin(samples_inc_2)]

onehot = OneHotEncoder(sparse=True).fit_transform(cell_subset_2[['sample']])
adata = adata[cell_subset_2.index]

pseudobulk_sum = adata.X.transpose() * onehot
pseudobulk_sum = pd.DataFrame.sparse.from_spmatrix(data=pseudobulk_sum, index=adata.var_names, columns=cell_subset_2['sample'].astype("category").cat.categories)

pseudobulk_sum.to_csv(pseudobulk_loc, sep="\t", index_label="gene")
