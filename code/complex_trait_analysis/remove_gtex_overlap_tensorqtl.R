library(vroom)
library(tidyverse)

celltype_eqtls <- vroom("results/static_eqtl_followup/qtl_sets/tensorqtl/original/signif_variant_gene_pairs.bed")
gwas_overlaps <- read.table("results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/tensorqtl_gwas_snplist.txt")$V1
overlap_eqtls <- vroom("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.full_gtex_overlap.bed",
                       col_names=c("CHR", "START", "STOP", "EB_ENSG", "EB_HGNC", 
                                   "EB_VARIANT_ID", "CELLTYPE", "CHR1", "START1", 
                                   "STOP1", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT"))

novel_celltype_eqtls_gwas <- celltype_eqtls %>%
  filter(EB_VARIANT_ID %in% gwas_overlaps) %>%
  anti_join(select(overlap_eqtls, c(EB_HGNC, EB_VARIANT_ID))) %>%
  select(EB_VARIANT_ID) %>%
  distinct() %>%
  write_tsv("results/cellregmap_eqtl_calling/eb_cellid/pseudobulk_tmm/basic/novel_tensorqtl_gwas_snplist.txt", col_names = F)
  