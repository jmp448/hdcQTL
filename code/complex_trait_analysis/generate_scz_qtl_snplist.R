library(vroom)
library(tidyverse)
library(data.table)

eb_bed_loc <- snakemake@input[['eb_bed']]
gtex_bed_loc <- snakemake@input[['gtex_bed']]
scz_sumstats_loc <- snakemake@input[['scz_sumstats']]

snplist_loc <- snakemake@output[['non_gtex_scz_eqtls']]
non_gtex_scz_eqtls_loc <- snakemake@output[['non_gtex_scz_eqtls_bed']]

eb_bed <- vroom(eb_bed_loc)
gtex_overlap_bed <- vroom(gtex_bed_loc,
                          col_names=c("CHR", "START", "STOP", "EB_ENSG", "EB_HGNC", 
                                      "EB_VARIANT_ID", "CELLTYPE", "CHR1", "START1", 
                                      "STOP1", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT"))

non_gtex_eqtls <- anti_join(eb_bed, gtex_overlap_bed, by=c("EB_VARIANT_ID", "EB_HGNC"))

scz_variants <- fread(scz_sumstats_loc)
scz_variants$V5 <- as.numeric(scz_variants$V5) # while not necessary for SCZ analysis, this script is used for other traits where fread reads the pvalue column as character
scz_hits <- scz_variants[V5 <= 5e-8, .(V1, V3)]

non_gtex_scz_eqtls <- non_gtex_eqtls %>%
  mutate(`#CHR`=str_replace(`#CHR`, "chr", "")) %>%
  inner_join(scz_hits, by=c("#CHR"="V1", "END"="V3")) %>%
  write_tsv(non_gtex_scz_eqtls_loc)

# save a snplist to check if any of these tag GTEx eQTLs
write(unique(non_gtex_scz_eqtls$EB_VARIANT_ID), snplist_loc)