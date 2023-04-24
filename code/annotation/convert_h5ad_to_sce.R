library(tidyverse)
library(zellkonverter)  # note that this package requires R version 4.1.0
library(SingleCellExperiment)

h5ad_dir <- snakemake@input[['h5ad']]
sce_dir <- snakemake@output[['sce']]

sce <- readH5AD(h5ad_dir, X_name='logcounts',
                use_hdf5=TRUE,
                layers=FALSE,
                uns=FALSE,
                obsm=FALSE,
                obsp=FALSE,
                verbose=TRUE)

saveRDS(sce, sce_dir)