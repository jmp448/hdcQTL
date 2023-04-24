library(Matrix)
library(Matrix.utils)
library(SingleCellExperiment)
library(tidyverse)
library(CelliD)

fca_loc <- snakemake@input[['fca']]
hvg_loc <- snakemake@input[['hvg']]
fca_embedding_loc <- snakemake@output[['fca_embedded']]
signatures_loc <- snakemake@output[['signatures']]

fca <- readRDS(fca_loc)
hvg_fca <- read_tsv(hvg_loc) %>% pull(gene)

fca <- RunMCA(fca, slot="logcounts", nmcs=50, features=hvg_fca)
fca_signatures <- GetGroupGeneSet(fca, group.by="celltype", n.features=100)

saveRDS(fca, fca_embedding_loc)
saveRDS(fca_signatures, signatures_loc)
