library(fastTopics)
library(tidyverse)

de_result_loc <- snakemake@input[["de_result"]]
context_loading_loc <- snakemake@output[['context_loading']]

load(de_result_loc)
L_mat <- as_tibble(fit_multinom$L) %>%
  write_tsv(context_loading_loc)