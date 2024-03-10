library(vroom)
library(tidyverse)
library(data.table)

dynamic_qtl_bed_loc <- snakemake@input[['dynamics_bed']]
gtex_bed_loc <- snakemake@input[['gtex_bed']]

snplist_loc <- snakemake@output[['non_gtex_dynamics']]
non_gtex_dynamics_loc <- snakemake@output[['non_gtex_dynamics_bed']]

dynamic_qtl_bed <- vroom(dynamic_qtl_bed_loc)
gtex_overlap_bed <- vroom(gtex_bed_loc,
                          col_names=c("CHR", "START", "STOP", "EB_ENSG", "EB_HGNC", 
                                      "EB_VARIANT_ID", "CELLTYPE", "CHR1", "START1", 
                                      "STOP1", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT"))

non_gtex_eqtls <- anti_join(dynamic_qtl_bed, gtex_overlap_bed, by=c("EB_VARIANT_ID", "EB_HGNC")) %>%
  write_tsv(non_gtex_dynamics_loc)

# save a snplist for the opentargets search
write(unique(non_gtex_eqtls$EB_VARIANT_ID), snplist_loc)