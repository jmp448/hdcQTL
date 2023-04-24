import numpy as np
import pandas as pd
import random
import scipy
import scanpy as sc
import scdrs

eb_loc = snakemake.input['eb']
eb_mca = snakemake.input['mca']
eb_celltype = snakemake.input['celltypes']

preprocessed_loc = snakemake.output['eb_preprocessed']

eb = sc.read_h5ad(eb_loc)

mca_embedding = pd.read_csv(eb_mca, sep="\t").set_index('cell').to_numpy()
eb.obsm['MCA'] = mca_embedding

# Load cell type labels & subset to cells successfully annotated
celltypes = pd.read_csv(eb_celltype, sep="\t").set_index('cell').reindex(eb.obs_names)
celltypes = celltypes[celltypes['value'] != 'unassigned']

### Subsample cell types with over 5K cells and remove cell types with fewer than 1K
celltypes = celltypes.groupby('value')
subsampled_labs = pd.DataFrame()
for name, group in celltypes:
    # If the number of observations in the group is between 1K and 5K,
    # add all of the observations to the subsampled dataframe
    if 1000 <= len(group) <= 5000:
        subsampled_labs = pd.concat([subsampled_labs, group])
    elif len(group) > 5000:
        # If the number of observations in the group is greater than the maximum,
        # randomly sample the maximum number of observations from the group
        subsampled_labs = pd.concat([subsampled_labs, group.sample(n=5000)])

eb = eb[subsampled_labs.index]
eb.obs['celltype'] = subsampled_labs['value']
eb.X = eb.layers['log1pPF']

scdrs.preprocess(eb, adj_prop='celltype')
sc.pp.neighbors(eb, use_rep='MCA')

eb.write_h5ad(preprocessed_loc)
