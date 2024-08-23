library(tidyverse)
library(vroom)

# snp_loc <- snakemake@input[[1]]
# bed_loc <- snakemake@output[[1]]

# snps <- vroom(snp_loc, col_select=c("feature_id", "snp_id")) %>%
#   rename(GTEX_GENE=feature_id) %>%
#   separate(snp_id, into=c("#CHR", "END", "GTEX_REF", "GTEX_ALT")) %>%
#   mutate(`#CHR`=paste0("chr", `#CHR`)) %>%
#   mutate(START=as.numeric(END)-1) %>%
#   select(c(`#CHR`, START, END, GTEX_GENE, GTEX_REF, GTEX_ALT)) %>%
#   write_tsv(bed_loc)

snp_locs <- snakemake@input
bed_loc <- snakemake@output[[1]]

snps <- vroom(snp_locs, col_select=c("chr", "snp_pos", "ref", "alt", "gene")) %>%
  rename(`#CHR`=chr, GTEX_GENE=gene, GTEX_REF=ref, GTEX_ALT=alt, END=snp_pos) %>%
  mutate(`#CHR`=paste0("chr", `#CHR`)) %>%
  mutate(START=as.numeric(END)-1) %>%
  select(c(`#CHR`, START, END, GTEX_GENE, GTEX_REF, GTEX_ALT)) %>%
  write_tsv(bed_loc)
