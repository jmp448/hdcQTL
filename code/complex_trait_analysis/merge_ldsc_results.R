library(vroom)
library(tissue)

ldsc_results_loc <- c("results/ldsc/baseline_eqtl_gtexpilotcomb_mash-signif_annot/UKB_460K.biochemistry_Testosterone_Male.results", "results/ldsc/baseline_eqtl_gtexpilotcomb_mash-signif_annot/UKB_460K.disease_HYPOTHYROIDISM_SELF_REP.results")

ldsc_results_loc <- snakemake@input

# helper functions
pull_trait <- function(s) {
  str_sub(tail(str_split(s, "/")[[1]], n=1), 1, -9)
}

# compile a list of QTL test p-values across all tests
results <- vroom(ldsc_results_loc, id="path") %>%
  mutate(trait=sapply(path, pull_trait), .keep="unused") %>%
  select()
  group_by(context) %>% 
  mutate(q=qvalue(pval_beta, lambda=0.85)$q) %>%
  write_tsv(all_qtl_loc)