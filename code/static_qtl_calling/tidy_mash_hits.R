library(mashr)
library(tidyverse)

# mash_model <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash_fitted_model.tophits.rds"

mash_model <- snakemake@input[['mash_model']]
hit_list <- snakemake@output[['mash_hits']]
lfsr_cutoff <- snakemake@wildcards[['threshold']]

m <- readRDS(mash_model)

sighits <- tibble(gv=names(get_significant_results(m, thresh=lfsr_cutoff))) %>%
  separate(gv, into=c("EB_HGNC", "EB_RSID"), sep="_") %>%
  write_tsv(hit_list)