library(tidyverse)
library(zellkonverter)  # note that this package requires R version 4.1.0
library(SingleCellExperiment)

eb_lowpass_h5ad_dir <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.TEMP.h5ad"
eb_lowpass_sce_dir <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.TEMP.sce"
# eb_lowpass_h5ad_dir <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.endoderm.raw.h5ad"
# eb_lowpass_sce_dir <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/Lowpass.3seqbatches.merged.endoderm.raw.sce"

eb_lowpass <- readH5AD(eb_lowpass_h5ad_dir)
# eb_lowpass <- readH5AD(eb_lowpass_h5ad_dir, X_name="counts")

# save the count data to the proper assay
eb_counts <- assay(eb_lowpass)
counts(eb_lowpass) <- t(eb_counts)
assay(eb_lowpass, "X") <- NULL

saveRDS(eb_lowpass, eb_lowpass_sce_dir)