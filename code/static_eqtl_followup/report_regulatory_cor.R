library(tidyverse)
library(vroom)

effects_loc <- snakemake@input[['effects']]
pos_loc <- snakemake@input[['pos']]

cor_loc <- snakemake@output[['cor']]

pos <- read_tsv(pos_loc, col_names="variant_id_rsid")
all_effects <- vroom(effects_loc) %>%
  filter(variant_id_rsid %in% pos$variant_id_rsid)

regulatory_cor <- cor(all_effects$EB_BHAT, all_effects$GTEX_BHAT)

write.table(regulatory_cor, file=cor_loc, quote=F, col.names=F, row.names=F)
