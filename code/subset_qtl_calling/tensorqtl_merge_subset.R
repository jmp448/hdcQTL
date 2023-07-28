library(tidyverse)
library(vroom)
library(qvalue)
set.seed(1234)

qtl_files <- snakemake@input
all_qtl_loc <- as.character(snakemake@output[['all_qtls']])
top_qtl_loc <- as.character(snakemake@output[['top_qtls']])

# helper functions
pull_tissue <- function(s) {
    str_extract(tail(str_split(s, "/")[[1]], n=1), "[^.]+")
}

# compile a list of QTL test p-values across all tests
qtls <- vroom(qtl_files, id="path") %>%
  mutate(context=sapply(path, pull_tissue), .keep="unused") %>%
  group_by(context) %>% 
  mutate(q=qvalue(pval_beta, lambda=0.85)$q) %>%
  write_tsv(all_qtl_loc)

sighits <- filter(qtls, q <= 0.05) %>%
  write_tsv(top_qtl_loc)
