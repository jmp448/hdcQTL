library(tidyverse)
library(qvalue)

eqtl_file <- snakemake@input[["eqtls"]]
all_qtls_file <- snakemake@output[["all_qtls"]]

eqtls <- read_tsv(eqtl_file) %>%
  add_count(EB_HGNC) %>%
  mutate(P_BONF=if_else(P_CELLREGMAP*n>1, 1, P_CELLREGMAP*n))

top_eqtls <- eqtls %>%
  group_by(EB_HGNC) %>%
  arrange(P_BONF) %>%
  slice_head(n=1)

top_eqtls$q <- qvalue(top_eqtls$P_BONF)$qvalues
global_cutoff <- filter(top_eqtls, q <= 0.1) %>%
  ungroup %>%
  arrange(desc(P_BONF)) %>%
  slice_head(n=1) %>% pull(P_BONF)

all_eqtls <- eqtls %>%
  filter(P_BONF <= global_cutoff)

write_tsv(all_eqtls, all_qtls_file)