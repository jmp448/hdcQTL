library(Matrix)
library(Matrix.utils)
library(scran)
library(SingleCellExperiment)
library(tidyverse)
library(CelliD)

# Read snakemake inputs
subsampled_loc <- snakemake@input[['fca_subsampled']]
hvg_loc <- snakemake@output[['hvgs']]
subsampled_norm_loc <- snakemake@output[['fca_subsampled_lognorm']]

# Load subsampled fetal cell atlas
fca_subsampled <- readRDS(subsampled_loc)

# Normalize data
libs <- librarySizeFactors(fca_subsampled)
sizeFactors(fca_subsampled) <- libs
fca_subsampled <- logNormCounts(fca_subsampled)

# Get highly variable features
fca_subsampled.dec <- modelGeneVar(fca_subsampled)
hvg.fca_subsampled <- getTopHVGs(fca_subsampled.dec, fdr.threshold=0.1)

tibble(gene=hvg.fca_subsampled) %>%
  write_tsv(hvg_loc)

saveRDS(fca_subsampled, subsampled_norm_loc)