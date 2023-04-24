import scanpy as sc
import pandas as pd
import numpy as np

scvi_embedded_adata_loc = snakemake.input[0]
umap_embedding_loc = snakemake.output[0]
umap_embedded_adata_loc = snakemake.output[1]

# Load data
adata = sc.read_h5ad(scvi_embedded_adata_loc)

# Get neighbors graph
sc.pp.neighbors(adata, use_rep="X_scVI")

# Get umap embedding
sc.tl.umap(adata)
np.savetxt(umap_embedding_loc, adata.obsm['X_umap'])

# Get clusters
sc.tl.leiden(adata, resolution=0.1)

# Save full object
adata.write_h5ad(umap_embedded_adata_loc)
