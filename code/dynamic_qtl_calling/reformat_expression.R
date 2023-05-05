###
# Write expression data to BED format for tensorQTL
###

library(tidyverse)

expression_loc <- snakemake@input[['exp']]
gene_loc <- snakemake@input[['gene_locs']]
expression_reformatted_loc <- snakemake@output[['exp']]

exp <- read_tsv(expression_loc)

genes <- read_tsv(gene_loc)

bed <- genes %>%
  rename(gene=hgnc, `#chr`=chr, start=`tss_start`, end=`tss_end`) %>%
  relocate(`#chr`, start, end, gene) %>%
  mutate(start=start-1) %>%
  inner_join(exp, by="gene") %>%
  write_tsv(expression_reformatted_loc)