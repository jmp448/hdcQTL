library(tidyverse)
library(vroom)

dynamic_egenes_loc <- snakemake@input[['egenes']]
all_tests_loc <- snakemake@input[['alltests']]
fdr_cutoff <- as.numeric(snakemake@wildcards[['fdr']])
sig_eqtl_loc <- snakemake@output[['sighits']]

dynamic_egenes <- vroom(dynamic_egenes_loc) %>%
  filter(pval_adj_bh <= fdr_cutoff) 
all_tests <- vroom(all_tests_loc)

# Get a genome-wide P-value cutoff (the Eigen-MT adjusted P-value of the gene closest to the FDR cutoff)
p_cutoff <- dynamic_egenes %>%
  slice_max(pval_adj_bh, n=1) %>%
  pull(pval_emt)

dynamic_eqtls_sig <- inner_join(all_tests, select(dynamic_egenes, c(phenotype_id, tests_emt)), by="phenotype_id") %>%
  mutate(pval_adj_emt=map2_dbl(pval_gi, tests_emt, `*`)) %>%
  filter(pval_adj_emt <= p_cutoff) %>%
  write_tsv(sig_eqtl_loc)