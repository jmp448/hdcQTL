library(Matrix)
library(Matrix.utils)
library(scran)
library(SingleCellExperiment)
library(tidyverse)
library(CelliD)
library(RcppArmadillo)

test_loc <- snakemake@input[['test_77']]
signatures_loc <- snakemake@input[['signatures_77']]

labeled_loc <- snakemake@output[['labeled_test_77']]

test_sce <- readRDS(test_loc)
fca_signatures <- readRDS(signatures_loc)

# Run hypergeometric test for cell type classification
cellid_enrichments <- RunCellHGT(test_sce, pathways = fca_signatures, dims = 1:50, n.features=100)

# Pull labels - highest enrichment w P<=0.01 in each cell
cellid_labels <- rownames(cellid_enrichments)[apply(cellid_enrichments, 2, which.max)]
cellid_labels <- ifelse(apply(cellid_enrichments, 2, max)>=2, yes = cellid_labels, "unassigned") 

# Save
test_sce$cellid_label <- cellid_labels
saveRDS(test_sce, labeled_loc)
