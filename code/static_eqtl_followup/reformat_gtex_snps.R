library(tidyverse)
library(vroom)

# snp_loc <- "/project2/gilad/jpopp/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.signif_variant_gene_pairs.txt"

snp_loc <- snakemake@input[[1]]
bed_loc <- snakemake@output[[1]]

snps <- vroom(snp_loc, col_select=c("variant_id", "gene_id")) %>%
  mutate(GTEX_GENE=str_extract(gene_id, "[^.]+")) %>%
  separate(variant_id, into=c("#CHR", "END", "GTEX_REF", "GTEX_ALT", "build")) %>%
  mutate(START=as.numeric(END)-1) %>%
  select(c(`#CHR`, START, END, GTEX_GENE, GTEX_REF, GTEX_ALT)) %>%
  write_tsv(bed_loc)
