library(tidyverse)
library(qusage)

plurinet_signatures_loc <- snakemake@input[['plurinet_signatures']]
fca_signatures_loc <- snakemake@input[['fca_signatures']]

combined_signatures_loc <- snakemake@output[['all_signatures']]

plurinet_gs <- read.delim(plurinet_signatures_loc)
plurinet_gs <- colnames(plurinet_gs)[3:ncol(plurinet_gs)]

fca_signatures <- readRDS(fca_signatures_loc)
fca_signatures$IPSC <- plurinet_gs

saveRDS(fca_signatures, combined_signatures_loc)
