import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import random

random.seed(1234)

adata_input_loc = snakemake.input['adata']
adata_output_loc = snakemake.output['adata']

adata = sc.read_h5ad(adata_input_loc)

# Remove batch
pilot_batches = ['P1','P2','P3']
adata_woPilot = adata[adata.obs['library_prep_batch'].isin(pilot_batches) == False, :]

# sampling cells from outlier individuals with too many cells represented by taking the median number of cells 
# across individuals
donors = set(adata_woPilot.obs.donor_id.to_list())
donor_cell = {}
for donor in donors:
    donor_cell[donor] = len(adata_woPilot.obs[adata_woPilot.obs['donor_id']==donor].index)

median_ncells = int(np.median(list(donor_cell.values())))
print("Median # cells per individual is {}".format(median_ncells))

# cell to donor mapping 
sample_mapping = adata_woPilot.obs[['donor_id']]

# individuals to be sampled: 
sample_donor = ["NA18499", "NA18508","NA18858","NA18862"]
# initialize the index from donors that are not sampled
new_donor_cell_index = sample_mapping.donor_id[sample_mapping.donor_id.isin(sample_donor)==False].index
for donor in sample_donor: 
    sample_index = sample_mapping.donor_id[sample_mapping.donor_id == donor].sample(n=median_ncells,replace=False,random_state=0).index
    # append sampled index to new_donor_cell_index
    new_donor_cell_index = new_donor_cell_index.append(sample_index)
    
# construct the new cell donor mapping
new_sample_mapping = sample_mapping.loc[new_donor_cell_index]

# subset the anndata object to the cells after sampling approach
adata_woPilot_sampled = adata_woPilot[new_sample_mapping.index, :]

adata_woPilot_sampled.write_h5ad(adata_output_loc)
