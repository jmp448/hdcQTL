library(tidyverse)
library(vroom)
library(qvalue) 
library(TreeBH)

ctprops_combined <- vroom("results/trans_qtl_calling/eb_cellid/pseudobulk_tmm/basic/all-celltype-proportions-sva/GOBP_TISSUE_DEVELOPMENT-variants.gtex_all_tissues.tsv") %>%
  mutate(tissue=str_extract(tissue, "[^.]+"))
variant_info <- vroom("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/trans_eqtl_variant_candidate_info.GOBP_TISSUE_DEVELOPMENT.Artery-Tibial.tsv")

# going to filter to the variants I grouped
mtc3 <- left_join(ctprops_combined, variant_info, by="variant_id") %>%
  group_by(tissue, phenotype_id, phenotype_id_ensg) %>%
  add_count(name="n_tests_per_locus") %>%
  mutate(locus.bonf.p=map2_dbl(pval, n_tests_per_locus, p.adjust, method="bonferroni"))

mtc2 <- mtc3 %>%
  slice_min(order_by=pval, n=1,  with_ties=F) %>%
  group_by(tissue, phenotype_id) %>%
  add_count(name="n_loci_per_celltype") %>%
  mutate(celltype.bonf.p = map2_dbl(locus.bonf.p, n_loci_per_celltype, p.adjust, method="bonferroni"))

mtc1 <- mtc2 %>%
  slice_min(order_by=pval, n=1, with_ties=F) %>%
  group_by(tissue) %>%
  mutate(tissue.bh = p.adjust(celltype.bonf.p, method="BH"))

celltypes_with_qtl <- filter(mtc1, tissue.bh <= 0.1) %>%
  mutate(rejected=T)
loci_with_qtl <- mtc2 %>% 
  left_join(select(celltypes_with_qtl, c(tissue, phenotype_id, rejected)), by=c("tissue", "phenotype_id"))

mtc1_tissue_level <- mtc2 %>%
  slice_min(order_by=pval, n=1, with_ties=F) %>%
  group_by(tissue) %>%
  mutate(tissue.bonf.p = p.adjust(celltype.bonf.p, method="bonferroni"))

mtc0_tissue_level <- mtc1_tissue_level %>%
  slice_min(order_by=pval, n=1, with_ties=F) %>%
  ungroup %>%
  mutate(global.bh = p.adjust(tissue.bonf.p, method="BH"))

# Compare power of TreeBH approach with just a bonferroni scan
power_comp <- left_join(ctprops_combined, variant_info, by="variant_id") %>%
  group_by(tissue) %>%
  add_count(name="n_tests_per_tissue") %>%
  mutate(alltests.bonf.p=map2_dbl(pval, n_tests_per_tissue, p.adjust, method="bonferroni"))
celltypes_with_qtl_notree <- filter(power_comp, alltests.bonf.p <= 0.1)

# Actually implement TreeBH approach
apply_treebh_in_tissue <- function(t, ctprops_all_tissues) {
  treebh_input_groups <- filter(ctprops_all_tissues, tissue==t) %>%
    group_by(phenotype_id) %>%
    mutate(celltype_id=cur_group_id()) %>%
    group_by(phenotype_id, phenotype_id_ensg) %>%
    mutate(egene_id=cur_group_id()) %>%
    group_by(phenotype_id, phenotype_id_ensg, variant_id) %>%
    mutate(test_id=cur_group_id()) 
  treebh_input_groups_df <- treebh_input_groups %>%
    ungroup() %>%
    select(celltype_id, egene_id, test_id) %>%
    as.data.frame()
  treebh_input_pvals <- treebh_input_groups %>%
    pull(pval)
  
  treebh_selections <- get_TreeBH_selections(treebh_input_pvals, treebh_input_groups_df, q=rep(0.5, 3), test="simes")
  
  treebh_selection <- 
}
treebh_input_full <- left_join(ctprops_combined, variant_info, by="variant_id")
treebh_input_groups <- treebh_input_full %>%
  group_by(tissue) %>%
  mutate(tissue_id = cur_group_id()) %>%
  group_by(tissue, phenotype_id) %>%
  mutate(celltype_id=cur_group_id()) %>%
  group_by(tissue, phenotype_id, phenotype_id_ensg) %>%
  mutate(egene_id=cur_group_id()) %>%
  group_by(tissue, phenotype_id, phenotype_id_ensg, variant_id) %>%
  mutate(test_id=cur_group_id()) 
treebh_input_groups_df <- treebh_input_groups %>%
  ungroup() %>%
  select(tissue_id, celltype_id, egene_id, test_id) %>%
  as.data.frame()
treebh_input_pvals <- treebh_input_groups %>%
  pull(pval)
um <- get_TreeBH_selections(treebh_input_pvals, treebh_input_groups_df, q=rep(0.25, 4), test="simes")

# See how many groups there are at each level
hist(treebh_input_full %>% group_by(tissue) %>% summarize(n=n_distinct(phenotype_id)) %>% pull(n))
hist(treebh_input_full %>% group_by(tissue, phenotype_id) %>% summarize(n=n_distinct(phenotype_id_ensg)) %>% pull(n))
hist(treebh_input_full %>% group_by(tissue, phenotype_id, phenotype_id_ensg) %>% summarize(n=n_distinct(variant_id)) %>% pull(n))

# Trans eQTL analysis
alltissues_combined <- vroom("results/trans_qtl_calling/eb_cellid/pseudobulk_tmm/basic/GOBP_CIRCULATORY_SYSTEM_PROCESS-genes_merged/GOBP_CIRCULATORY_SYSTEM_DEVELOPMENT-variants.gtex_all_tissues.tsv") %>%
  mutate(tissue=str_extract(tissue, "[^.]+"))
circsystem_tissues <- c("Artery_Coronary", "Artery_Aorta", "Artery_Tibial", "Heart_Atrial_Appendage", "Heart_Left_Ventricle")
circsystem_combined <- filter(alltissues_combined, tissue %in% circsystem_tissues)

mtc1_trans <- left_join(circsystem_combined, variant_info, by="variant_id") %>%
  group_by(tissue, phenotype_id, phenotype_id_ensg) %>%
  add_count(name="n_tests_per_locus") %>%
  mutate(locus.bonf.p=map2_dbl(pval, n_tests_per_locus, p.adjust, method="bonferroni"))
mtc2_trans_end <- mtc1_trans %>%
  slice_min(order_by=pval, n=1,  with_ties=F) %>%
  group_by(tissue, phenotype_id) %>%
  mutate(celltype.bh = p.adjust(locus.bonf.p, method="BH"))
trans_sighits_tissue_level <- filter(mtc2_trans_end, celltype.bh <= 0.1)

common_trans_sighits <- trans_sighits_tissue_level %>%
  group_by(phenotype_id_ensg, phenotype_id) %>%
  add_count() %>%
  filter(n >= 2) # focus on the trans effects detected in multiple tissues

mtc2_trans <- mtc1_trans %>%
  slice_min(order_by=pval, n=1,  with_ties=F) %>%
  group_by(tissue, phenotype_id) %>%
  add_count(name="n_loci_per_celltype") %>%
  mutate(celltype.bonf.p = map2_dbl(locus.bonf.p, n_loci_per_celltype, p.adjust, method="bonferroni"))

mtc3_trans <- mtc2_trans %>%
  slice_min(order_by=pval, n=1, with_ties=F) %>%
  group_by(tissue) %>%
  mutate(tissue.bh = p.adjust(celltype.bonf.p, method="BH"))

trans_sighits <- filter(mtc3_trans, tissue.bh <= 0.1)


