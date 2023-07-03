library(tidyverse)
library(vroom)
library(mashr)

# mash_loc="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash_fitted_model.full.rds"
# harmonized_tests_loc="results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt"
# bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/Retinal-cells/genotypes_filtered_plink.bim"
# gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
#bed_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.bed"

mash_loc <- snakemake@input[['mash_model']]
harmonized_tests_loc <- snakemake@input[['harmonized_tests']]
# bim_file <- snakemake@input[['bim_file']]
# gtf_loc <- snakemake@input[['gtf_loc']]

bed_file <- snakemake@output[['bedfile']]

# Harmonized tests
harmonized_tests <- vroom(harmonized_tests_loc) 

# Get mash significant hits
m <- readRDS(mash_loc)
mash_hits <- as_tibble(get_lfsr(m)[get_significant_results(m, thresh=0.05),], rownames="test") %>%
  pivot_longer(!test, names_to="EB_CELLTYPE", values_to="lfsr") %>%
  filter(lfsr <= 0.05) %>%
  group_by(test) %>%
  summarize(EB_CELLTYPE=paste(EB_CELLTYPE, collapse=",")) %>%
  separate(test, into=c("EB_HGNC", "EB_VARIANT_ID"), sep="_") %>%
  left_join(harmonized_tests, by=c("EB_VARIANT_ID"="variant_id_rsid", "EB_HGNC"="phenotype_id_hgnc")) %>%
  drop_na() %>% # this will remove tests not run in GTEx
  dplyr::rename(EB_ENSG=phenotype_id_ensg) %>%
  separate(variant_id, into=c("#CHR", "END", "EB_REF", "EB_ALT", "BUILD")) %>%
  mutate(END=as.numeric(END), START=END-1) %>%
  select(c(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID, EB_CELLTYPE)) %>%
  write_tsv(bed_file)