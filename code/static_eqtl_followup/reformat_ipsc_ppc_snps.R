library(tidyverse)
library(vroom)
library(stats)

snp_locs <- snakemake@input[[1]]
bed_loc <- snakemake@output[[1]]

load(ppc_eqtls)

# subset to significant eQTLs 
hits <- as.data.table(gene)[, .SD[which.min(bonferroni)], by=gene_id]
hits$bh = p.adjust(hits$bonferroni)
snps <- gene %>%
  rename(`#CHR`=chrom, GTEX_GENE=gene_id, GTEX_REF=ref, GTEX_ALT=alt, END=pos) %>%
  mutate(GTEX_GENE=str_extract(GTEX_GENE, "[^.]+")) %>%
  mutate(`#CHR`=paste0("chr", `#CHR`)) %>%
  mutate(START=as.numeric(END)-1) %>%
  select(c(`#CHR`, START, END, GTEX_GENE, GTEX_REF, GTEX_ALT)) %>%
  write_tsv(bed_loc)
