library(tidyverse)
library(vroom)

all_gtex_loc <- snakemake@input[['all_gtex']]
all_eb_tests_loc <- snakemake@input[['all_eb']]
gtf_loc <- snakemake@input[['gtf']]
rsid_snpinfo_map_loc <- snakemake@input[['rsid_map']]
harmonized_tests_loc <- snakemake@output[['harmonized_tests']]

all_gtex_tests <- vroom(all_gtex_loc, col_names=c("phenotype_id_ensg", "variant_id"))
all_eb_tests <- vroom(all_eb_tests_loc) %>%
  rename(phenotype_id_hgnc=phenotype_id, variant_id_rsid=variant_id) %>%
  select(c(variant_id_rsid, phenotype_id_hgnc)) %>%
  distinct()
rsid_snpinfo_map <- vroom(rsid_snpinfo_map_loc, col_names=c("variant_id_rsid", "variant_id"))

# Subset SNP info map to tests run in the EBs, with unique rsIDs
rsid_snpinfo_map <- rsid_snpinfo_map %>%
  filter(variant_id_rsid %in% unique(all_eb_tests$variant_id_rsid)) %>%
  filter(variant_id_rsid != ".")

# Reformat EB variant IDs
all_eb_tests <- all_eb_tests %>%
  left_join(rsid_snpinfo_map, by=c("variant_id_rsid")) 

hgnc_ensg_dict <- vroom(gtf_loc, col_select=c("hgnc", "ensg"))

all_eb_tests <- all_eb_tests %>%
  left_join(hgnc_ensg_dict, by=c("phenotype_id_hgnc"="hgnc")) %>%
  drop_na(ensg) %>% dplyr::rename(phenotype_id_ensg=ensg)

# Find overlap with GTEx, accounting for the possibility of allele flips
asis_overlap <- inner_join(all_eb_tests, all_gtex_tests, by=c("variant_id", "phenotype_id_ensg"))

# double check there are no allele flips
eb_variants_alleleflipped <- all_eb_tests %>%
  separate(variant_id, into=c("CHR", "POS", "REF", "ALT", "BUILD"), sep="_") %>%
  unite("variant_id", c(CHR, POS, ALT, REF, BUILD), sep="_") %>%
  pull(variant_id)
stopifnot(length(intersect(eb_variants_alleleflipped, all_gtex_tests$variant_id)) == 0)

# save the dictionary of common (between EBs and the tissues selected) tests
tests_harmonized <- write_tsv(asis_overlap, harmonized_tests_loc)