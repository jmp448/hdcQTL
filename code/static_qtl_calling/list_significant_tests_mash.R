library(tidyverse)
library(vroom)
library(mashr)

# permutation_results_loc <- "results/static_qtl_calling/eb_cmstages/pseudobulk_tmm/basic/8pcs/tensorqtl_permutations.all.tsv"
# nominal_loc <- "results/static_qtl_calling/eb_cmstages/pseudobulk_tmm/basic/8pcs/tensorqtl_nominal.all.tsv"

mash_loc <- snakemake@input[['full_fitted_model']]
significant_variant_gene_pairs_loc <- snakemake@output[['hit_list']]
  
m <- readRDS(mash_loc)
mash_hits <- as_tibble(get_lfsr(m)[get_significant_results(m, thresh=0.05),], rownames="test") %>%
  pivot_longer(!test, names_to="EB_CELLTYPE", values_to="lfsr") %>%
  filter(lfsr <= 0.05) %>%
  group_by(test) %>%
  summarize(EB_CELLTYPE=paste(EB_CELLTYPE, collapse=",")) %>%
  separate(test, into=c("EB_HGNC", "EB_VARIANT_ID"), sep="_") %>%
  write_tsv(significant_variant_gene_pairs_loc)