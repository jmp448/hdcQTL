import scanpy as sc
import pandas as pd

full_adata_loc = snakemake.input[0]
fca_hvg_loc = snakemake.input[1]
subset_adata_loc = snakemake.output[0]

# Load EB data
adata = sc.read_h5ad(full_adata_loc)

# Load FCA highly variable gene set
fca_hvgs = pd.read_csv(fca_hvg_loc, sep="\t")
overlap_hvgs = fca_hvgs.loc[fca_hvgs['gene'].isin(adata.var_names)]

# Filter EB data to FCA hvgs and save
adata_hvg = adata[:, overlap_hvgs['gene']]
adata_hvg.X = adata_hvg.layers['log1pPF']

adata_hvg.write_h5ad(subset_adata_loc)