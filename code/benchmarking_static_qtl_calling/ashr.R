library(tidyverse)
library(ashr)
library(vroom)
set.seed(1234)

eqtls_loc <- snakemake@input[['eqtls']]
random_tests_loc <- snakemake@input[['random']]
top_hits_loc <- snakemake@input[['top']]
df_loc <- snakemake@input[['df']]

trained_model_loc <- snakemake@output[['trained']]
tophit_fitted_loc <- snakemake@output[['tophits']]

# compute standard errors, omitting tests where we can't obtain reasonable estimates
qtls <- vroom(eqtls_loc, col_select=c(SNP, gene, beta, `t-stat`, df)) %>%
  mutate(se=if_else(as.logical((beta==0) * (`t-stat`==0)), rep(NA_real_, nrow(.)), beta/`t-stat`), .after='beta') %>%
  drop_na() %>%
  unite(gv, c(gene, SNP), sep="--", remove=FALSE)
df.shared <- read.table(df_loc)$V1

# learn the prior, g, on a subset of tests
random_tests <- read_tsv(random_tests_loc, col_names="test") %>% pull(test)
random_qtls <- qtls %>%
  filter(gv %in% random_tests)

betahat_random <- random_qtls %>%
  select(gv, beta) %>%
  deframe

sehat_random <- random_qtls %>%
  select(gv, se) %>%
  deframe

ash_trained <- ash(betahat_random, sehat_random, df=df.shared)
saveRDS(ash_trained, trained_model_loc)

# test on the strongest test per gene
top_hits <- read_tsv(top_hits_loc, col_names="test") %>% pull(test)
top_qtls <- qtls %>%
  filter(gv %in% top_hits)

betahat_tophits <- top_qtls %>%
  select(gv, beta) %>%
  deframe

sehat_tophits <- top_qtls %>%
  select(gv, se) %>%
  deframe

ash_test <- ash(betahat_tophits, sehat_tophits, df=df.shared, g=get_fitted_g(ash_trained), fixg=T)
saveRDS(ash_test, tophit_fitted_loc)