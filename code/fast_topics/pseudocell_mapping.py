import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import random

random.seed(1234)

adata_input_loc = snakemake.input['adata_subsampled']
cell_pseudocell_mapping_loc = snakemake.output['pseudocell_mapping']

# adapted from Ben Strober
def assign_pseudocells_inplace(adata, cluster_resolution, batch_label):
    adata.obs['donor_batch']= adata.obs['donor_id'].astype(str)+'_'+adata.obs[batch_label].astype(str)
    unique_donor_batches = sorted(np.unique(adata.obs['donor_batch']))
    num_cells = len(adata.obs['donor_batch'])
    # Initialize vector to keep track of cluster assignments
    cluster_assignments = np.asarray(['unassigned']*num_cells,dtype='<U40')

    # Loop through individuals
    for db in unique_donor_batches:
        # Get cell indices corresponding to this individual
        cell_indices = adata.obs['donor_batch'] == db
        # Number of cells in this indiviudal
        num_cells_per_indi_batch = sum(cell_indices)
        # Make anndata object for just this individual
        adata_db = adata[cell_indices, :]
        # Construct neighborhood graph for cells from this indivudal
        sc.pp.neighbors(adata_db, use_rep="X_scVI")
        # Perform leiden clustering
        sc.tl.leiden(adata_db, resolution=cluster_resolution)
        # Get leiden cluster assignemnts
        leiden_cluster_assignments = adata_db.obs['leiden']
        # Add to global vector of assignments
        cluster_assignments[cell_indices] = np.char.add(db + '_', leiden_cluster_assignments.astype(str))
        # Delete adata_indi from memory
        del adata_db
    adata.obs['pseudocell_' + str(cluster_resolution)] = cluster_assignments

adata = sc.read_h5ad(adata_input_loc)
assign_pseudocells_inplace(adata, cluster_resolution=15, batch_label='collection.date')

cluster_assignments = adata.obs[["pseudocell_15"]]
cluster_assignments.to_csv(cell_pseudocell_mapping_loc, sep='\t')


