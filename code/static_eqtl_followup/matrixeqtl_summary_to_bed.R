library(tidyverse)

eqtl_file <- snakemake@input[[1]]
bed_file <- snakemake@output[[1]]

#eqtl_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.tsv"
#bed_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.bed"

bed <- vroom(eqtl_file) %>%
  filter(q <= 0.1) %>%
  select(c(SNP, gene)) %>%
  dplyr::rename(GENE = gene) %>%
  mutate(`#CHR`=str_extract(SNP, "[^_]+")) %>%
  mutate(END=as.integer(str_extract(SNP, "[^_]+$")), .keep="unused") %>%
  mutate(START=as.integer(END - 1)) %>%
  relocate(`#CHR`, START, END, GENE) %>%
  write_tsv(bed_file)
  