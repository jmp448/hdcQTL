library(tidyverse)
library(vroom)

eqtl_file <- snakemake@input[['qtl_summary']]
bim_file <- snakemake@input[['bim_file']]
pos_file <- snakemake@output[['posfile']]
celltype <- as.character(snakemake@wildcards[['type']])

#eqtl_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.tsv"
#bed_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.bed"

bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS", "variant_id"))

pos <- vroom(eqtl_file) %>%
  filter(context==celltype) %>%
  dplyr::rename(GENE=phenotype_id) %>%
  select(c(variant_id, GENE)) %>%
  left_join(bim, by="variant_id") %>%
  select(-c(variant_id)) %>%
  relocate(`#CHR`, POS, GENE) %>%
  write_tsv(pos_file)
  