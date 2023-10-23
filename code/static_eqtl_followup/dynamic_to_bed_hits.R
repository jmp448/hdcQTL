library(tidyverse)
library(vroom)
library(mashr)

qtl_loc <- snakemake@input[['eqtls']]
bim_file <- snakemake@input[['bim_file']]
gtf_loc <- snakemake@input[['gtf_loc']]
traj <- snakemake@wildcards[['trajectory']]

bed_file <- snakemake@output[['bedfile']]

## Get a mapping from rsID to chromosome positions
### Note that some variants will be missing from this file bc they only pass MAF threshold 
### after subsetting donors in one of the celltypes that didn't contain everyone
bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS_BP", "variant_id")) %>%
  dplyr::rename(END=`POS_BP`) %>%
  mutate(START=as.integer(END)-1) # bed files are zero-indexed

# Get a mapping from HGNC to ENSG
gencode <- vroom(gtf_loc) %>% 
  select(c(hgnc, ensg))

# Merge tests from all three trajectories
qtls <- vroom(qtl_loc)
qtls_bed <- qtls %>%
  filter(variant_id != ".") %>%
  dplyr::rename(GENE=phenotype_id) %>%
  select(c(variant_id, GENE)) %>%
  distinct() %>%
  inner_join(bim, by="variant_id") %>%
  inner_join(gencode, by=c("GENE"="hgnc")) %>%
  dplyr::rename(EB_ENSG=ensg, EB_HGNC=GENE, EB_VARIANT_ID=variant_id) %>%
  relocate(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID) %>%
  mutate(EB_TRAJECTORY=traj) %>%
  write_tsv(bed_file)
