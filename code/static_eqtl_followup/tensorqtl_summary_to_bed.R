library(tidyverse)
library(vroom)

eqtl_file <- snakemake@input[['qtl_summary']]
bim_file <- snakemake@input[['bim_file']]
bed_file <- snakemake@output[['bedfile']]
celltype <- as.character(snakemake@wildcards[['type']])

#eqtl_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.tsv"
#bed_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.bed"

bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS_BP", "variant_id")) %>%
  dplyr::rename(END=`POS_BP`) %>%
  mutate(START=as.integer(END)-1) # bed files are zero-indexed

bed <- vroom(eqtl_file) %>%
  filter(context==celltype) %>%
  dplyr::rename(GENE=phenotype_id) %>%
  select(c(variant_id, GENE)) %>%
  left_join(bim, by="variant_id") %>%
  select(-c(variant_id)) %>%
  relocate(`#CHR`, START, END, GENE) %>%
  write_tsv(bed_file)
  