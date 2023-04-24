import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import scdrs

adata_loc = snakemake.input['adata']
preprocessed_anndata_loc = snakemake.output['adata_preprocessed']

# Load and preprocess EB anndata
eb_adata = sc.read_h5ad(adata_loc)
eb_adata.X = eb_adata.layers['log1pPF']
del eb_adata.layers['counts']
del eb_adata.layers['PFlog1pPF']
scdrs.preprocess(eb_adata)
eb_adata.write_h5ad(preprocessed_anndata_loc)
