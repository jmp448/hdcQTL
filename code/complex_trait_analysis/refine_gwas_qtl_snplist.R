library(vroom)
library(tidyverse)

eb_bed_loc <- snakemake@input[['non_gtex_gwas_eqtls']]
gtex_bed_loc <- snakemake@input[['gtex_overlap_tag_variants']]

refined_bed_loc <- snakemake@output[['filtered_bed']]

non_gtex_gwas_eqtls_bed <- vroom(eb_bed_loc)

tagged_gtex_bed <- vroom(gtex_bed_loc,
                          col_names=c("CHR", "START", "STOP", "EB_ENSG", "EB_HGNC", 
                                      "EB_VARIANT_ID", "CELLTYPE", "CHR1", "START1", 
                                      "STOP1", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT"))
non_gtex_gwas_eqtls_notags_bed <- anti_join(non_gtex_gwas_eqtls_bed, tagged_gtex_bed, by=c("EB_HGNC", "EB_VARIANT_ID")) %>%
  write_tsv(refined_bed_loc)


