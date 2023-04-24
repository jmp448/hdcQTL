import scanpy as sc
import pandas as pd
import numpy as np

raw_counts_loc = snakemake.input[0]
adata_qc_loc = snakemake.output[0]

adata = sc.read_h5ad(raw_counts_loc)
sc.pp.calculate_qc_metrics(adata, inplace=True)
adata.write_h5ad(adata_qc_loc)
