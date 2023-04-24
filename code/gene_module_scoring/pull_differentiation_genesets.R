library(tidyverse)

gs_full <- read_tsv("/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.bp.Hs.symbols.gs")

is_dev_module <- function(trait) {
  grepl("DIFFERENTIATION", trait, fixed=T) + grepl("DEVELOPMENT", trait, fixed=T) > 0
}

gs_full <- gs_full %>%
  rowwise() %>%
  mutate(is_dev=is_dev_module(TRAIT)) %>%
  mutate(n_genes=str_count(GENESET, ",") + 1)

gs_modsize <- gs_full %>%
  filter((n_genes >= 50) & (n_genes <= 200)) 

gs_differentiation <- gs_modsize %>%
  filter(is_dev == T) %>%
  select(TRAIT, GENESET) %>%
  write_tsv("/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.bp.Hs.symbols.differentiation.gs")