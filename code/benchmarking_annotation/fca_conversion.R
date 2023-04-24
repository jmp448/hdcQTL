library(tidyverse)
library(zellkonverter)

fca_full <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/counts.subsampled.mca_pca_fromhvgs.sce")

writeH5AD(fca_full, "/project2/gilad/jpopp/ebQTL/data/fca/counts.subsampled.mca_pca_fromhvgs.h5ad",
          X_name='logcounts')