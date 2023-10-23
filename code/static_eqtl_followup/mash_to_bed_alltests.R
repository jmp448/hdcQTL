library(tidyverse)
library(vroom)
library(mashr)

# mash_loc="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash_fitted_model.full.rds"
# bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/PNS-Glia/genotypes_filtered_plink.bim"
# gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
#bed_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.bed"

mash_loc <- snakemake@input[['mash_model']]
bim_file <- snakemake@input[['bim_file']]
gtf_loc <- snakemake@input[['gtf_loc']]

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

# Get mash significant hits
m <- readRDS(mash_loc)
mash_tests <- tibble(test=rownames(get_lfsr(m))) %>%
  separate(test, into=c("EB_HGNC", "EB_VARIANT_ID"), sep="_") %>%
  filter(EB_VARIANT_ID != ".") %>%
  inner_join(bim, by=c("EB_VARIANT_ID"="variant_id")) %>%
  inner_join(gencode, by=c("EB_HGNC"="hgnc")) %>%
  dplyr::rename(EB_ENSG=ensg) %>%
  select(c(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID)) %>%
  write_tsv(bed_file)