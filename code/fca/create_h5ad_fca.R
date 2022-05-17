library(tidyverse)
library(Matrix)
library(reticulate)
use_condaenv("scvi-gpflow")
use_python("/home/jpopp/.conda/envs/scvi-gpflow/bin/python", required=TRUE)
library(anndata)

# Save the inputs from snakemake
counts_file <- snakemake@input[[1]]
cell_file <- snakemake@input[[2]]
gene_file <- snakemake@input[[3]]
output_h5ad_file <- snakemake@output[[1]]
output_mapping_file <- snakemake@output[[2]]

# For debugging
counts_file <- "data/fca/counts.sampled.rds"
cell_file <- "data/fca/cell_metadata.rds"
gene_file <- "data/fca/gene_metadata.rds"
output_h5ad_file <- "data/single_cell_objects/fca.sampled.h5ad"
output_mapping_file <- "data/fca/organ_celltype.tsv"
  
# Load data - counts, cell & gene metadata
counts <- readRDS(counts_file) %>%
  t

gene_metadata <- readRDS(gene_file) %>%
  filter(gene_id %in% colnames(counts))

cell_metadata <- readRDS(cell_file) %>%
  filter(sample %in% rownames(counts))

# Separate the UMAP embedding from the cell metadata
umap <- cell_metadata %>%
  select(Main_cluster_umap_1, Main_cluster_umap_2) %>%
  as.matrix

cell_metadata <- cell_metadata %>%
  select(-c(Main_cluster_umap_1, Main_cluster_umap_2))

# Create anndata object
fca_adata <- AnnData(
  X=counts,
  obs=cell_metadata,
  var=gene_metadata,
  obsm=list(X_umap=umap)
)

write_h5ad(fca_adata, output_file)

# Additionally, it will be useful to create a mapping of cell types to organ of origin
organ_celltype <- select(cell_metadata, Organ_cell_lineage) %>%
  unique %>%
  separate(Organ_cell_lineage, into=c("organ", "celltype"), sep="-") %>%
  write_tsv(output_mapping_file)
