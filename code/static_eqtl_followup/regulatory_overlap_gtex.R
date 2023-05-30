library(tidyverse)
library(vroom)

# all_gtex_loc <- "results/static_eqtl_followup/gtex/allpairs_filtered/Liver.allpairs.filtered_50kb_dist2tss.txt"
# hits_gtex_loc <- "results/static_eqtl_followup/gtex/sigtests_filtered/Liver.signif_variant_gene_pairs.filtered_50kb_dist2tss.txt"
# all_eb_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/tensorqtl_nominal.all.tsv"
# hits_eb_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.tsv"
# harmonized_tests_loc <- "results/static_eqtl_followup/gtex/allpairs_filtered/eb_gtex_harmonized_tests.txt"

all_gtex_loc <- snakemake@input[['all_gtex']]
hits_gtex_loc <- snakemake@input[['hits_gtex']]
all_eb_loc <- snakemake@input[['all_eb']]
hits_eb_loc <- snakemake@input[['hits_eb']]
harmonized_tests_loc <- snakemake@input[['harmonized_tests']]

tissue <- as.character(snakemake@wildcards[['tissue']])
type <- as.character(snakemake@wildcards[['celltype']])

positions_loc <- snakemake@output[['joint_hit_positions']]
effects_loc <- snakemake@output[['joint_effects']]
  
harmonized_tests <- vroom(harmonized_tests_loc)
harmonized_variants <- harmonized_tests %>%
  select(variant_id, variant_id_rsid) %>%
  distinct()

# Load GTEx data and filter to tests that were shared between studies
all_gtex <- vroom(all_gtex_loc, col_select = c("gene_id", "variant_id", "slope", "pval_nominal")) %>%
  rename(phenotype_id_ensg=gene_id) %>% 
  inner_join(select(harmonized_tests, c(variant_id, phenotype_id_ensg)), by=c("variant_id", "phenotype_id_ensg"))
hits_gtex <- vroom(hits_gtex_loc, col_select = c("gene_id", "variant_id", "slope")) %>%
  rename(phenotype_id_ensg=gene_id) %>% 
  inner_join(select(harmonized_tests, c(variant_id, phenotype_id_ensg)), by=c("variant_id", "phenotype_id_ensg"))

all_eb <- vroom(all_eb_loc) %>%
  filter(celltype==type) %>%
  rename(phenotype_id_hgnc=phenotype_id, variant_id_rsid=variant_id) %>% 
  inner_join(harmonized_tests, by=c("variant_id_rsid", "phenotype_id_hgnc")) %>%
  select(-c(variant_id_rsid, phenotype_id_hgnc, celltype))
hits_eb <- vroom(hits_eb_loc) %>%
  filter(celltype==type) %>%
  rename(phenotype_id_hgnc=phenotype_id, variant_id_rsid=variant_id) %>% 
  inner_join(harmonized_tests, by=c("variant_id_rsid", "phenotype_id_hgnc")) %>%
  select(-c(variant_id_rsid, phenotype_id_hgnc, celltype))

# Get the union of all hits
all_hits <- select(hits_gtex, c(variant_id, phenotype_id_ensg)) %>%
  bind_rows(select(hits_eb, c(variant_id, phenotype_id_ensg))) %>%
  distinct()

# Save the positions in one file, as rs IDs, for clumping
pos <- all_hits %>%
  select(variant_id) %>%
  distinct() %>%
  left_join(harmonized_variants, by="variant_id") %>%
  select(variant_id_rsid) %>%
  write_tsv(positions_loc, col_names=F)

# Save the estimated effect sizes from each study at all positions to another file
eb_effects <- all_eb %>%
  inner_join(all_hits, by=c("variant_id", "phenotype_id_ensg")) %>%
  select(variant_id, phenotype_id_ensg, slope) %>%
  dplyr::rename(EB_BHAT=slope) %>%
  drop_na()

gtex_effects <- all_gtex %>%
  inner_join(all_hits, by=c("variant_id", "phenotype_id_ensg")) %>%
  select(variant_id, phenotype_id_ensg, slope) %>%
  dplyr::rename(GTEX_BHAT=slope) %>%
  drop_na()

all_effects <- inner_join(eb_effects, gtex_effects, by=c("variant_id", "phenotype_id_ensg")) %>%
  left_join(harmonized_variants, by="variant_id") %>%
  write_tsv(effects_loc)