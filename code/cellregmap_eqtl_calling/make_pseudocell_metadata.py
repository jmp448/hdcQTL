import scanpy as sc
import numpy as np
import pandas as pd
import xarray as xr

pseudocell_mapping_loc = snakemake.input['pseudocell_map']
cell_metadata_loc = snakemake.input['metadata']

pseudocell_metadata_loc = snakemake.output['pseudocell_metadata'] 

## import the pseudocell anndata object (raw counts)
metadata = pd.read_csv(cell_metadata_loc)[['Line.True', 'sex']].drop_duplicates().rename(columns={'Line.True': 'donor_id'})
metadata['donor_id'] = ['NA' + str(s) for s in metadata['donor_id']]
metadata['sex'] = np.where(metadata['sex'] == 'F', 0, 1)

cell_map = pd.read_csv(pseudocell_mapping_loc, sep="\t").rename(columns={'individual_leiden_clusters_15': 'pseudocell'})
cell_map['donor_id'] = cell_map['pseudocell'].str.split('_').str[0]
pseudocell_map = cell_map[['donor_id', 'pseudocell']].drop_duplicates()

pseudocell_map = pseudocell_map.merge(metadata, on='donor_id', how='inner').sort_values(by=["donor_id", "pseudocell"])

pseudocell_map.to_csv(pseudocell_metadata_loc, sep = "\t")
