library(tidyverse)
library(vroom)
library(qvalue) 
library(TreeBH)

ctprops_combined <- vroom("results/trans_qtl_calling/eb_cellid/pseudobulk_tmm/basic/all-celltype-proportions/GOBP_TISSUE_DEVELOPMENT-variants.gtex_all_tissues.tsv") %>%
  mutate(tissue=str_extract(tissue, "[^.]+")) %>%
  rename(celltype=phenotype_id, variant_id_gtex=variant_id) 
variant_info <- vroom("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/trans_eqtl_variant_candidate_info.GOBP_TISSUE_DEVELOPMENT.tsv") %>%
  rename(egene=phenotype_id)

mtc <- left_join(ctprops_combined, variant_info, by="variant_id_gtex") %>%
  group_by(tissue, celltype) %>%
  slice_min(order_by=pval_beta, n=1, with_ties=F) %>%
  group_by(tissue) %>%
  mutate(tissue.bh = p.adjust(pval_perm, method="BH"))

ctprop_qtls <- filter(mtc, tissue.bh <= 0.1)
ctprop_maybes <- filter(ctprops_combined, pval_beta <= 0.05)
# effective_tests <- vroom("results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/trans_eqtl_variant_effective_tests.GOBP_TISSUE_DEVELOPMENT.tsv") %>%
#   dplyr::rename(gene=`...1`)
# 
# # going to filter to the variants I grouped
# mtc_gene_level <- left_join(ctprops_combined, variant_info, by="variant_id_gtex") %>%
#   left_join(effective_tests, by=c("egene"="gene")) %>%
#   drop_na() %>% # figure out why there's an NA
#   group_by(tissue, celltype, egene) %>%
#   mutate(locus.bonf.p=map2_dbl(pval, n_eff, p.adjust, method="bonf"))
# 
# mtc_celltype_level <- mtc_gene_level %>%
#   slice_min(order_by=pval, n=1,  with_ties=F) %>%
#   group_by(tissue, celltype) %>%
#   add_count(name="n_loci_per_celltype") %>%
#   mutate(celltype.bonf.p = map2_dbl(locus.bonf.p, n_loci_per_celltype, p.adjust, method="holm"))
# 
# mtc_tissue_level <- mtc_celltype_level %>%
#   slice_min(order_by=pval, n=1, with_ties=F) %>%
#   group_by(tissue) %>%
#   mutate(tissue.bh = p.adjust(celltype.bonf.p, method="BH"))
# 
# celltypes_with_qtl <- filter(mtc_tissue_level, tissue.bh <= 0.1) %>%
#   mutate(rejected=T)
# 
# # Compare power of TreeBH approach with just a bonferroni scan for each celltype
# mtc_celltype_level_2 <- left_join(ctprops_combined, variant_info, by="variant_id_gtex") %>%
#   group_by(tissue, celltype) %>%
#   mutate(n_eff = 77) %>%
#   mutate(p.adj_test=map2_dbl(pval, n_eff, p.adjust, method="bonf"))
# 
# mtc_tissue_level_2 <- mtc_celltype_level_2 %>%
#   slice_min(order_by=p.adj_test, n=1, with_ties=F) %>%
#   group_by(tissue) %>%
#   mutate(p.adj_celltype = p.adjust(p.adj_test, method="BH"))
# 
# celltypes_with_qtl <- filter(mtc_tissue_level, tissue.bh <= 0.1) %>%
#   mutate(rejected=T)
# 
# # And compare to an EigenMT-informed scan
# apply_eigenmt_control <- function(p, m_eff) {
#   return(min(p * m_eff, 1))
# }
# 
# mtc_celltype_level_3 <- left_join(ctprops_combined, variant_info, by="variant_id_gtex") %>%
#   group_by(tissue, celltype) %>%
#   add_count(name="n_tests_per_celltype") %>%
#   mutate(n_eff_per_celltype=99) %>%
#   mutate(p.adj_test= map2_dbl(pval, n_eff_per_celltype, apply_eigenmt_control))
# 
# mtc_tissue_level_3 <- mtc_celltype_level_3 %>%
#   slice_min(order_by=p.adj_test, n=1, with_ties=F) %>%
#   group_by(tissue) %>%
#   mutate(p.adj_celltype = p.adjust(p.adj_test, method="BH"))
# 
# celltypes_with_qtl <- filter(mtc_tissue_level, tissue.bh <= 0.1) %>%
#   mutate(rejected=T)
# 
# 
# # # Actually implement TreeBH approach
# # apply_treebh_in_tissue <- function(t, ctprops_all_tissues) {
# #   treebh_input_groups <- filter(ctprops_all_tissues, tissue==t) %>%
# #     group_by(phenotype_id) %>%
# #     mutate(celltype_id=cur_group_id()) %>%
# #     group_by(phenotype_id, phenotype_id_ensg) %>%
# #     mutate(egene_id=cur_group_id()) %>%
# #     group_by(phenotype_id, phenotype_id_ensg, variant_id) %>%
# #     mutate(test_id=cur_group_id()) 
# #   treebh_input_groups_df <- treebh_input_groups %>%
# #     ungroup() %>%
# #     select(celltype_id, egene_id, test_id) %>%
# #     as.data.frame()
# #   treebh_input_pvals <- treebh_input_groups %>%
# #     pull(pval)
# #   
# #   treebh_selections <- get_TreeBH_selections(treebh_input_pvals, treebh_input_groups_df, q=rep(0.5, 3), test="simes")
# #   
# #   treebh_selection <- 
# # }
# # treebh_input_full <- left_join(ctprops_combined, variant_info, by="variant_id")
# # treebh_input_groups <- treebh_input_full %>%
# #   group_by(tissue) %>%
# #   mutate(tissue_id = cur_group_id()) %>%
# #   group_by(tissue, phenotype_id) %>%
# #   mutate(celltype_id=cur_group_id()) %>%
# #   group_by(tissue, phenotype_id, phenotype_id_ensg) %>%
# #   mutate(egene_id=cur_group_id()) %>%
# #   group_by(tissue, phenotype_id, phenotype_id_ensg, variant_id) %>%
# #   mutate(test_id=cur_group_id()) 
# # treebh_input_groups_df <- treebh_input_groups %>%
# #   ungroup() %>%
# #   select(tissue_id, celltype_id, egene_id, test_id) %>%
# #   as.data.frame()
# # treebh_input_pvals <- treebh_input_groups %>%
# #   pull(pval)
# # um <- get_TreeBH_selections(treebh_input_pvals, treebh_input_groups_df, q=rep(0.1, 4), test="simes")
# # 
# # # See how many groups there are at each level
# # hist(treebh_input_full %>% group_by(tissue) %>% summarize(n=n_distinct(phenotype_id)) %>% pull(n))
# # hist(treebh_input_full %>% group_by(tissue, phenotype_id) %>% summarize(n=n_distinct(phenotype_id_ensg)) %>% pull(n))
# # hist(treebh_input_full %>% group_by(tissue, phenotype_id, phenotype_id_ensg) %>% summarize(n=n_distinct(variant_id)) %>% pull(n))
# 
# # what if we just used the top QTL per gene (as opposed to applying EigenMT)
# eb_signif <- vroom("/project2/gilad/jpopp/ebQTL/results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.tsv")
# eb_top_any_celltype <- eb_signif %>%
#   group_by(phenotype_id) %>%
#   slice_min(order_by=pval_nominal, n=1, with_ties=F) 
# variants_top <- inner_join(variant_info, dplyr::select(eb_top_any_celltype, c(variant_id, phenotype_id)), by=c("variant_id"="variant_id", "egene"="phenotype_id"))
# 
# mtc_tests <- left_join(ctprops_combined, variants_top, by="variant_id_gtex") %>%
#   drop_na() %>%
#   group_by(tissue, celltype) %>%
#   add_count(name="n_tests_per_celltype") %>%
#   mutate(p.adj_tests=map2_dbl(pval, n_tests_per_celltype, p.adjust, method="hommel"))
# 
# mtc_celltypes <- mtc_tests %>%
#   slice_min(order_by=p.adj_tests, n=1, with_ties=F) %>%
#   group_by(tissue) %>%
#   mutate(p.adj_celltypes = p.adjust(p.adj_tests, method="BH"))
# 
# # to get a hint for how much EigenMT or a permutation scheme could help, skip the bonferroni step across tests
# mtc_alltests_pseudoeigen <- left_join(ctprops_combined, variant_info, by="variant_id_gtex") %>%
#   group_by(tissue, celltype, egene) %>%
#   add_count(name="n_tests_per_egene") %>%
#   mutate(adj.p.tests=33*pval)
#   mutate(adj.p.tests=map2_dbl(pval, n_tests_per_egene, p.adjust, method="none"))
# 
# mtc_genes_pseudoeigen <- mtc_alltests_pseudoeigen %>%
#   slice_min(order_by=pval, n=1, with_ties=F) %>%
#   group_by(tissue, celltype) %>%
#   add_count(name="n_genes_per_celltype") %>%
#   mutate(adj.p.genes=map2_dbl(pval, n_genes_per_celltype, p.adjust, method="hommel"))
# 
# mtc_celltypes_pseudoeigen <- mtc_genes_pseudoeigen %>%
#   slice_min(order_by=pval, n=1, with_ties=F) %>%
#   group_by(tissue) %>%
#   mutate(tissue.bh = p.adjust(adj.p.genes, method="BH"))
# 
# # More targeted testing - nervous system development gene set
# brain_tissues <- c("Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", 
#                  "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9",
#                  "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia",
#                  "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra") 
# ctprops_brain <- filter(ctprops_combined, tissue %in% brain_tissues)
# gmt_file <- "/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
# 
# gmt_lines <- readLines(gmt_file)
# 
# gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, "\t")))
# gmt_list2 <- lapply(gmt_list, function(x) list(x[1], x[2], x[3:length(x)]))
# 
# gmt_df <- as_tibble(do.call("rbind", gmt_list2)) %>%
#   mutate(across(V1:V2, unlist)) %>%
#   dplyr::rename(geneset=V1, link=V2, genes=V3) %>%
#   filter(str_sub(geneset, 1, 4) == "GOBP") # filter to BP 
# cns_genes <- unlist(filter(gmt_df, geneset == "GOBP_CENTRAL_NERVOUS_SYSTEM_DEVELOPMENT")$genes)
# cns_qtls <- filter(variant_info, egene %in% cns_genes)
# ctprops_brain_cns <- filter(ctprops_brain, variant_id_gtex %in% cns_qtls$variant_id_gtex)
# 
# mtc_brain <- left_join(ctprops_brain_cns, cns_qtls, by="variant_id_gtex") %>%
#   group_by(tissue, celltype) %>%
#   add_count(name="n_genes_per_celltype") %>%
#   mutate(locus.bonf.p=map2_dbl(pval, n_genes_per_celltype, p.adjust, method="bonf"))
# 
# ## Muscle
# muscle_tissues <- c("Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Muscle_Skeletal") 
# ctprops_muscle <- filter(ctprops_combined, tissue %in% muscle_tissues)
# 
# muscle_genes <- unlist(filter(gmt_df, geneset == "GOBP_MUSCLE_TISSUE_DEVELOPMENT")$genes)
# muscle_qtls <- filter(variant_info, egene %in% muscle_genes)
# ctprops_muscle <- filter(ctprops_muscle, variant_id_gtex %in% muscle_qtls$variant_id_gtex)
# 
# mtc_muscle <- left_join(ctprops_muscle, muscle_qtls, by="variant_id_gtex") %>%
#   group_by(tissue, celltype) %>%
#   add_count(name="n_genes_per_celltype") %>%
#   mutate(locus.bonf.p=map2_dbl(pval, n_genes_per_celltype, p.adjust, method="bonf"))