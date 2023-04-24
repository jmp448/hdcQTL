library(tidyverse)

eqtl_file <- snakemake@input[[1]]
positions_file <- snakemake@output[[1]]

#eqtl_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.tsv"
#positions_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.pos"

pos <- vroom(eqtl_file) %>%
  filter(q <= 0.1) %>%
  select(c(SNP, gene)) %>%
  dplyr::rename(GENE = gene) %>%
  mutate(`#CHR`=str_extract(SNP, "[^_]+")) %>%
  mutate(POS=as.integer(str_extract(SNP, "[^_]+$")), .keep="unused") %>%
  relocate(`#CHR`, POS, GENE) %>%
  write_tsv(positions_file)