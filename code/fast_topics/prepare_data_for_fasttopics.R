#! /usr/bin/env Rscript

# Install package 'anndata' to read anndata (h5ad)
# Install package 'sceasy' to convert h5ad to seurat object (.rds)
# Required modules: R, gcc - 'module load R gcc'
# Installation: devtools::install_github("cellgeni/sceasy"); ref: https://github.com/cellgeni/sceasy 

# Load a few packages.
library(Matrix)
library(Seurat)
library(stringi)

# Load the dataset
full <- readRDS(snakemake@input[[1]])
# Get the raw counts data for the full dataset
counts<- GetAssayData(full, slot="counts", assay="RNA")
counts<- t(counts)

cat(sprintf("Number of samples:%d\n", nrow(counts)))
cat(sprintf("Number of genes:%d\n", ncol(counts)))
cat(sprintf("Proportion of counts that are non-zero: %0.1f%%.\n",
            100*mean(counts > 0)))

# Save the counts and metadata

save(list = c("genes", "counts"),
     file= snakemake@output[[1]])

