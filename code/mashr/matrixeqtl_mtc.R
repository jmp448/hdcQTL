library(tidyverse)
library(qvalue)

eqtl_file <- snakemake@input[["eqtls"]]
df_file <- snakemake@input[["df"]]
all_tests_file <- snakemake@output[["all_tests"]]
top_hits_file <- snakemake@output[["top_tests"]]
n_hits_file <- snakemake@output[["n_hits"]]

df.obs <- read_tsv(df_file, col_names=F)$X1[[1]]

eqtls <- read_tsv(eqtl_file) %>%
  mutate(df=df.obs, .after=`t-stat`) %>%
  left_join(count(., gene), by="gene") %>%
  mutate(bonf.p=if_else(`p-value`*n>1, 1, `p-value`*n)) %>%
  select(-n) %>%
  write_tsv(all_tests_file)

top_eqtls <- eqtls %>%
  group_by(gene) %>%
  arrange(bonf.p) %>%
  filter(!duplicated(gene))

top_eqtls$q <- qvalue(top_eqtls$bonf.p)$qvalues

write_tsv(top_eqtls, top_hits_file)

write(sum(top_eqtls$q <= 0.1), file=n_hits_file)
