library(tidyverse)
library(vroom)
library(qvalue)
set.seed(1234)

# Get input & output file names from snakemake
# cell.types <- read_tsv("data/benchmark_specificity_methods/eb_cellid/eb_cellid.pseudobulk_tmm.tsv", n_max=0, col_select=-c(gene)) %>%
#   colnames %>%
#   sapply(function(s){str_sub(s, start=7)}) %>%
#   unique()
# qtl_files <- paste0("results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/basic/", cell.types, "/8pcs/matrixeqtl.cis_qtl_pairs.all.tsv")

qtl_files <- snakemake@input
all_qtl_loc <- as.character(snakemake@output[['all_qtls']])
top_qtl_loc <- as.character(snakemake@output[['top_qtls']])

# helper functions
pull_type <- function(s) {
  str_split(s, "/")[[1]][6]
}

# compile a list of QTL test p-values across all tests
qtls <- vroom(qtl_files, id="path") %>%
  mutate(context=sapply(path, pull_type), .keep="unused") %>%
  group_by(context) %>% 
  mutate(q=qvalue(pval_beta, lambda=0.85)$q) %>%
  write_tsv(all_qtl_loc)

sighits <- filter(qtls, q <= 0.05) %>%
  write_tsv(top_qtl_loc)
