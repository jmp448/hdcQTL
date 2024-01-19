library(tidyverse)
library(qusage)

plurinet_gs <- read.delim("data/gene_sets/MUELLER_PLURINET.v2023.2.Hs.gmt")
plurinet_gs <- colnames(plurinet_gs)[3:ncol(plurinet_gs)]

fca_signatures <- readRDS("data/fca/fca_signatures.rds")
fca_signatures$IPSC <- plurinet_gs

saveRDS(fca_signatures, "data/fca/fca_signatures.with_ipsc.rds")
