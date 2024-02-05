library(Matrix)
library(Matrix.utils)
library(scran)
library(SingleCellExperiment)
library(tidyverse)
library(CelliD)

set.seed(42)

fca_counts_loc="data/fca/counts.subsampled.sce"

# Read snakemake inputs
sce_loc <- snakemake@input[['sce']]

train_subset_loc <- snakemake@output[['train_subset']]
test_subset_loc <- snakemake@output[['test_subset']]
signatures_subset_loc <- snakemake@output[['signatures_subset']]

# Load subsampled fetal cell atlas data
fca <- readRDS(fca_counts_loc)

# Subset half the data for training, leaving half for testing
train_cells <- sample(colnames(fca), size=round(ncol(fca)/2), replace=F)
test_cells <- setdiff(colnames(fca), train_cells)

fca_train <- fca[,train_cells]
fca_test <- fca[,test_cells]

# Normalize
libs <- librarySizeFactors(fca_train)
sizeFactors(fca_train) <- libs
fca_train <- logNormCounts(fca_train)

libs_test <- librarySizeFactors(fca_test)
sizeFactors(fca_test) <- libs_test
fca_test <- logNormCounts(fca_test)

fca_train.dec <- modelGeneVar(fca_train)
hvg.fca_train <- getTopHVGs(fca_train.dec, fdr.threshold=0.1)

# Learn signatures
fca_train <- RunMCA(fca_train, slot="logcounts", nmcs=50, features=hvg.fca_train)
fca_test <- RunMCA(fca_test, slot="logcounts", nmcs=50, features=hvg.fca_train)
fca_subset_signatures <- GetGroupGeneSet(fca_train, group.by="celltype", n.features=100)

saveRDS(fca_train, train_subset_loc)
saveRDS(fca_test, test_subset_loc)
saveRDS(fca_subset_signatures, signatures_subset_loc)

