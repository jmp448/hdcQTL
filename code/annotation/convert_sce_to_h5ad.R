library(tidyverse)
library(zellkonverter)  # note that this package requires R version 4.1.0
library(SingleCellExperiment)

sce_dir <- snakemake@output[['sce']]
h5ad_dir <- snakemake@input[['h5ad']]

sce <- readRDS(sce_dir)

writeH5AD(sce, h5ad_dir, 
          X_name='logcounts',
          verbose=TRUE)