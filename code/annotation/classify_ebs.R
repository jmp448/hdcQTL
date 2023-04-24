library(Matrix)
library(Matrix.utils)
library(SingleCellExperiment)
library(tidyverse)
library(CelliD)

eb_loc <- snakemake@input[['eb_hvgs']]
signatures_loc <- snakemake@input[['signatures']]

mca_embedding_loc <- snakemake@output[['mca_embedding']]
cellid_labels_loc <- snakemake@output[['cellid_labels']]
cellid_all_enrichments_loc <- snakemake@output[['cellid_all_enrichments']]

ebs <- readRDS(eb_loc)
fca_signatures <- readRDS(signatures_loc)

# Get an MCA embedding of the EB data
ebs <- RunMCA(ebs, slot="logcounts", nmcs=50)

# Run hypergeometric test for cell type classification
cellid_enrichments <- RunCellHGT(ebs, pathways = fca_signatures, dims = 1:50, n.features=100)

# Pull labels - highest enrichment w P<=0.01 in each cell
cellid_labels <- rownames(cellid_enrichments)[apply(cellid_enrichments, 2, which.max)]
cellid_labels <- ifelse(apply(cellid_enrichments, 2, max)>=2, yes = cellid_labels, "unassigned") 

# Save MCA embedding to a tsv file
mca_embedding <- as_tibble(reducedDim(ebs, "MCA"), rownames="cell") %>%
  write_tsv(mca_embedding_loc)

# Write all assignment information to a tsv file
cellid_enrichments %>%
  as.matrix %>% t %>%
  as_tibble(rownames="cell") %>%
  write_tsv(cellid_all_enrichments_loc)

# Write labels to a tsv file
as_tibble(cellid_labels, rownames="cell") %>%
  write_tsv(cellid_labels_loc)
