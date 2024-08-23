import os
import pandas as pd
from anndata import read_h5ad
import scdrs
import scanpy as sc

# read input files
scdrs_adata_loc = snakemake.input['adata_scdrs']
connectivities_adata_loc = snakemake.input['adata_conn']

# read output files
output_loc = snakemake.output['scdrs_conn']

# load data
scdrs = read_h5ad(scdrs_adata_loc)
conn = read_h5ad(connectivities_adata_loc)

# copy over cell connectivity graph
scdrs.obsp['connectivities']=conn.obsp['connectivities']

# save
scdrs.write_h5ad(output_loc)