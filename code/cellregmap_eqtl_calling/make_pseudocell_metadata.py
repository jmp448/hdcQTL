import scanpy as sc
import numpy as np
import pandas as pd
import xarray as xr

pseudocell_mapping_loc = snakemake.input['pseudocell_map']
cell_metadata_loc = snakemake.input['metadata']

pseudocell_metadata_loc = snakemake.output['pseudocell_metadata'] 

## import the pseudocell anndata object (raw counts)
metadata = pd.read_csv(cell_metadata_loc)[['Line.True', 'collection.date', 'sex']].drop_duplicates().rename(columns={'Line.True': 'donor_id', 'collection.date': 'batch'})
metadata = metadata.dropna() 
metadata['donor_id'] = ['NA' + str(s) for s in metadata['donor_id']]
metadata['sex'] = np.where(metadata['sex'] == 'F', 0, 1)
metadata['batch'] = [b.replace("/", ".") for b in list(metadata['batch'])]
batch_encoded = pd.get_dummies(metadata['batch'], prefix='batch')
metadata = pd.concat([metadata, batch_encoded], axis=1)

cell_map = pd.read_csv(pseudocell_mapping_loc, sep="\t").rename(columns={'pseudocell_15': 'pseudocell'})
cell_map['donor_id'] = cell_map['pseudocell'].str.split('_').str[0]
cell_map['batch'] = cell_map['pseudocell'].str.split('_').str[1]
cell_map['batch'] = [b.replace("/", ".") for b in list(cell_map['batch'])]

pseudocell_map = cell_map[['pseudocell', 'donor_id', 'batch']].drop_duplicates()
pseudocell_map = pseudocell_map.merge(metadata, on=['donor_id', 'batch'], how='inner').sort_values(by=["donor_id", "batch", "pseudocell"]).drop(columns="batch")

pseudocell_map.to_csv(pseudocell_metadata_loc, sep = "\t", index=False)
