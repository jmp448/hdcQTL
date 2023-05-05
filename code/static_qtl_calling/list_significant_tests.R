library(tidyverse)
library(vroom)

# permutation_results_loc <- "results/static_qtl_calling/eb_cmstages/pseudobulk_tmm/basic/8pcs/tensorqtl_permutations.all.tsv"
# nominal_loc <- "results/static_qtl_calling/eb_cmstages/pseudobulk_tmm/basic/8pcs/tensorqtl_nominal.all.tsv"

permutation_results_loc <- snakemake@input[['permutations']]
nominal_loc <- snakemake@input[['nominal']]
significant_variant_gene_pairs_loc <- snakemake@output[['hit_list']]
  
permutations <- vroom(permutation_results_loc)
nominal <- vroom(nominal_loc)

# First, identify the empirical p-value of the gene closest to the q-value cutoff
per_context_cutoffs <- permutations %>%
  group_by(context) %>%
  arrange(desc(q)) %>%
  filter(q<=0.05) %>%
  slice_head(n=1) %>%
  select(pval_perm, context) %>%
  rename(global_p_cutoff=pval_perm)

# Next, calculate a nominal P-value threshold for each gene based on the beta distribution parameters
per_gene_cutoffs <- permutations %>%
  left_join(per_context_cutoffs, by="context") %>%
  mutate(gene_cutoff = pmap_dbl(list(global_p_cutoff, beta_shape1, beta_shape2), qbeta)) %>%
  select(phenotype_id, context, gene_cutoff)

# Finally, list all variant-gene pairs which fall below this threshold
significant_tests <- nominal %>%
  left_join(per_gene_cutoffs, by=c("celltype"="context", "phenotype_id"="phenotype_id")) %>%
  filter(pval_nominal <= gene_cutoff) %>%
  write_tsv(significant_variant_gene_pairs_loc)