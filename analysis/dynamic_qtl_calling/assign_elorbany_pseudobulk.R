library(tidyverse)

# Save the pseudobulk 
elorbany_pseudobulk <- read_tsv("/project2/gilad/jpopp/sc-dynamic-eqtl/data/pseudobulk-cm/bin16/counts.tsv") %>%
  write_tsv("/project2/gilad/jpopp/ebQTL/data/dynamic_qtl_calling/elorbany-cm_16bins/pseudobulk_tmm/elorbany-cm_16bins.pseudobulk_tmm.tsv")

# Save the metadata
sample_counts <- read_tsv("/project2/gilad/jpopp/sc-dynamic-eqtl/results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin16/bin_ncells.tsv")
sample_pseudotimes <- read_tsv("/project2/gilad/jpopp/sc-dynamic-eqtl/results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin16/bin_medians.tsv")
sample_depths <- read_tsv("/project2/gilad/jpopp/sc-dynamic-eqtl/results/eqtl_dynamic/linear_dQTL/pseudobulk-cm/bin16/bin_libsize.tsv")

sample_summary <- sample_counts %>%
  dplyr::rename(ind_type=binind, n_cells_unfiltered=n) %>%
  mutate(individual=str_extract(ind_type, "[^_]+"), type=str_extract(ind_type, "[^_]+$"), n_cells_filtered=n_cells_unfiltered) %>%
  left_join(sample_pseudotimes, by=join_by(ind_type==binind)) %>%
  dplyr::rename(pseudotime=t) %>%
  left_join(select(sample_depths, c(sample, lib.size)), by=join_by(ind_type==sample)) %>%
  dplyr::rename(total_counts=lib.size) %>%
  mutate(dropped=is.na(total_counts)) %>%
  replace_na(list(total_counts=0)) %>%
  relocate(ind_type, individual, type, n_cells_unfiltered, dropped, total_counts, n_cells_filtered, pseudotime) %>%
  write_tsv("/project2/gilad/jpopp/ebQTL/data/dynamic_qtl_calling/elorbany-cm_16bins/pseudobulk_tmm/sample_summary.tsv")
