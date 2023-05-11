import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import random

random.seed(1234)

adata_input_loc = snakemake.input['adata_subsampled']
cell_pseudocell_mapping_loc = snakemake.output['pseudocell_mapping']

# adapted from Ben Strober
def perform_leiden_clustering_in_each_individual(adata, cluster_resolution):
    unique_individuals = sorted(np.unique(adata.obs['donor_id']))
    num_cells = len(adata.obs['donor_id'])
    # Initialize vector to keep track of cluster assignments
    cluster_assignments = np.asarray(['unassigned']*num_cells,dtype='<U40')

    # Loop through individuals
    for individual in unique_individuals:
        # Get cell indices corresponding to this individual
        cell_indices = adata.obs['donor_id'] == individual
        # Number of cells in this indiviudal
        num_cells_per_indi = sum(cell_indices)
        # Make anndata object for just this individual
        adata_indi = adata[cell_indices, :]
        # Construct neighborhood graph for cells from this indivudal
        sc.pp.neighbors(adata_indi, use_rep="X_scVI")
        # Perform leiden clustering
        sc.tl.leiden(adata_indi, resolution=cluster_resolution)
        # Get leiden cluster assignemnts
        leiden_cluster_assignments = adata_indi.obs['leiden']
        # Add to global vector of assignments
        cluster_assignments[cell_indices] = np.char.add(individual + '_', leiden_cluster_assignments.astype(str))
        # Delete adata_indi from memory
        del adata_indi
    adata.obs['individual_leiden_clusters_' + str(cluster_resolution)] = cluster_assignments
    return adata

adata = sc.read_h5ad(adata_input_loc)
perform_leiden_clustering_in_each_individual(adata, cluster_resolution=15)
cluster_assignments = adata.obs[["individual_leiden_clusters_15"]]
cluster_assignments.to_csv(cell_pseudocell_mapping_loc, sep='\t')
