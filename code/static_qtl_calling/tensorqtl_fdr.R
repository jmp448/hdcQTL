library(tidyverse)
library(vroom)
library(qvalue)

cis_df_loc <- as.character(snakemake@input)
cis_df_fdr_loc <- as.character(snakemake@output)

cis_df_fdr <- vroom(cis_df_loc) %>%
  mutate(qval=qvalue(pval_beta, lambda=0.85)$q) %>%
  write_tsv(cis_df_fdr_loc)
