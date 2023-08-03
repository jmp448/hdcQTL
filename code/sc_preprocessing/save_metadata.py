import scanpy as sc
import pandas as pd
import numpy as np

raw_counts_loc = snakemake.input[0]
metadata_loc = snakemake.output[0]

adata = sc.read_h5ad(raw_counts_loc)
adata.obs.to_csv(metadata_loc, sep="\t")
