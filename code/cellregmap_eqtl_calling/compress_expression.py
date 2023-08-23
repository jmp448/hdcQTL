import scanpy as sc
import numpy as np
import pandas as pd
import xarray as xr

pseudocell_adata_loc = snakemake.input['pseudocell_adata']
compressed_expression_loc = snakemake.output['exp']
tabular_expression_loc = snakemake.output['tabular_exp']

## import the pseudocell anndata object (raw counts)
adata = sc.read_h5ad(pseudocell_adata_loc)

# normalize pseudocell expression
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.normalize_total(adata)

# make the sparse matrix dense
adata.X = adata.X.toarray()

# make phenotype file as dataframe
phenotype = pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names)
phenotype.to_csv(tabular_expression_loc, sep="\t")

#Save phenotype as binary file for fast read-in
phenotype = xr.DataArray(
    data=phenotype.values,
    dims=["pseudo_cell", "gene"],
    coords={"pseudo_cell": phenotype.index.values, "gene": phenotype.columns.values},
    name="expression"
    )
# assert all(phenotype.pseudo_cell.values == sample_mapping.index.values)

phenotype.to_netcdf(compressed_expression_loc)
