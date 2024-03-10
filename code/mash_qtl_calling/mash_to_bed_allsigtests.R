library(tidyverse)
library(vroom)
library(mashr)

# mash_loc="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash_fitted_model.full.rds"
# harmonized_tests_loc="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt"
# bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/Retinal-cells/genotypes_filtered_plink.bim"
# gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
#bed_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.bed"

mash_loc <- snakemake@input[['mash_model']]
bim_file <- snakemake@input[['bim_file']]
gtf_loc <- snakemake@input[['gtf_loc']]

bed_file <- snakemake@output[['bedfile']]

# Get variant info
bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS_BP", "variant_id")) %>%
  dplyr::rename(END=`POS_BP`) %>%
  mutate(START=as.integer(END)-1) # bed files are zero-indexed

# Get a mapping from HGNC to ENSG
gencode <- vroom(gtf_loc) %>% 
  select(c(hgnc, ensg)) %>%
  dplyr::rename(EB_ENSG=ensg)

# Get mash significant hits
m <- readRDS(mash_loc)
mash_hits <- as_tibble(get_lfsr(m)[get_significant_results(m, thresh=0.05),], rownames="test") %>%
  pivot_longer(!test, names_to="EB_CELLTYPE", values_to="lfsr") %>%
  filter(lfsr <= 0.05) %>%
  group_by(test) %>%
  summarize(EB_CELLTYPE=paste(EB_CELLTYPE, collapse=",")) %>%
  separate(test, into=c("EB_HGNC", "EB_VARIANT_ID"), sep="_") %>%
  filter(EB_VARIANT_ID != ".") %>%
  inner_join(bim, by=c("EB_VARIANT_ID"="variant_id")) %>%
  inner_join(gencode, by=c("EB_HGNC"="hgnc")) %>%
  select(c(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID, EB_CELLTYPE)) %>%
  write_tsv(bed_file)