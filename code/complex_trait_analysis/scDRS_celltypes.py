import os
import pandas as pd
from anndata import read_h5ad
import scdrs
import scanpy as sc

# read input files
adata_loc = snakemake.input['adata']
celltype_label_loc = snakemake.input['celltypes']
scores_loc = snakemake.input['scores']

# read output files
output_loc = snakemake.output['celltype_scores']

# load anndata object, cell type labels, and scores
adata = read_h5ad(adata_loc)
celltype_labels = pd.read_csv(celltype_label_loc, sep="\t").set_index("cell")
assert adata.obs.index.equals(celltype_labels.index), "Cell type indices do not match anndata object"
adata.obs['celltype'] = celltype_labels['value']

scores = pd.read_csv(scores_loc, sep="\t").set_index("cell")
assert adata.obs.index.equals(scores.index), "Disease score indices do not match anndata object"

# run scDRS cell type analysis
celltype_stats = scdrs.method.downstream_group_analysis(
    adata=adata,
    df_full_score=scores,
    group_cols=["celltype"],
)["celltype"]

# save output
celltype_stats.to_csv(output_loc, sep="\t")
