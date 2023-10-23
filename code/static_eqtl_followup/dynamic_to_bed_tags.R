library(tidyverse)
library(vroom)
library(mashr)

qtl_loc <- snakemake@input[['eqtls']]
tags_loc <- snakemake@input[['tags']]
bim_file <- snakemake@input[['bim_file']]
traj <- snakemake@wildcards[['trajectory']]

bed_file <- snakemake@output[['bedfile']]

## Get a mapping from rsID to chromosome positions
### Note that some variants will be missing from this file bc they only pass MAF threshold 
### after subsetting donors in one of the celltypes that didn't contain everyone
bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS_BP", "variant_id")) %>%
  dplyr::rename(END=`POS_BP`) %>%
  mutate(START=as.integer(END)-1) # bed files are zero-indexed

# Create a sorted BED file 
qtls <- vroom(qtl_loc, col_select=c("EB_ENSG", "EB_HGNC", "EB_VARIANT_ID"))
tags_bed <- read.table(tags_loc, header=T) %>% 
  as_tibble() %>%
  separate_rows(TAGS) %>%
  select(TAGS, SNP) %>%
  right_join(qtls, by=c("SNP"="EB_VARIANT_ID"), multiple="all") %>% # get QTL info from the tagged variant
  rename(EB_TAGGED_SNP=SNP, EB_VARIANT_ID=TAGS) %>%
  inner_join(bim, by=c("EB_VARIANT_ID"="variant_id")) %>%
  relocate(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID, EB_TAGGED_SNP) %>%
  write_tsv(bed_file)
