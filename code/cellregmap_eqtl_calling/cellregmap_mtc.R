library(tidyverse)
library(qvalue)

eqtl_file <- snakemake@input[["eqtls"]]
all_tests_file <- snakemake@output[["all_tests"]]
top_hits_file <- snakemake@output[["top_tests"]]
sig_hits_file <- snakemake@output[["sig_hits"]]

eqtls <- read_tsv(eqtl_file) %>%
  add_count(EB_HGNC) %>%
  mutate(P_BONF=if_else(P_CELLREGMAP*n>1, 1, P_CELLREGMAP*n)) %>%
  select(-n) %>%
  write_tsv(all_tests_file)

top_eqtls <- eqtls %>%
  group_by(EB_HGNC) %>%
  arrange(P_BONF) %>%
  slice_head(n=1)

top_eqtls$q <- qvalue(top_eqtls$P_BONF)$qvalues
write_tsv(top_eqtls, top_hits_file)

sig_eqtls <- filter(top_eqtls, q <= 0.1) %>%
  write_tsv(sig_hits_file)