library(vroom)
library(qvalue)
library(tidyverse)

sarkar_sumstats_loc <- snakemake@input[[1]]
sarkar_sumstats_posfile <- snakemake@output[[1]]

sarkar_sumstats <- vroom(sarkar_sumstats_loc)
sarkar_sumstats$q <- qvalue(sarkar_sumstats$p_beta, lambda=0.85)$q

sarkar_hits <- sarkar_sumstats %>%
  filter(q <= 0.05) %>%
  select(chr, end, gene) %>%
  mutate(start=as.numeric(end)-1, .after=1) %>%
  rename(c("#CHR"="chr", "START"="start", "END"="end", "GTEX_GENE"="gene")) %>%
  mutate(GTEX_REF=NA_character_, GTEX_ALT=NA_character_) %>%
  write_tsv(sarkar_sumstats_posfile)