library(tidyverse)
library(mashr)
library(flashier)
set.seed(1234)

beta_df_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/tensorqtl_nominal.betas.tsv"
se_df_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/tensorqtl_nominal.standard_errors.tsv"
sample_summary_loc <- "data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/sample_summary_manual.tsv"

beta_df_loc <- snakemake@input[['beta_df']]
se_df_loc <- snakemake@input[['se_df']]
sample_summary_loc <- snakemake@input[['sample_summary']]

npcs <- as.numeric(snakemake@wildcards[['npcs']])
n_contexts_cutoff <- as.numeric(snakemake@params[['min_contexts']])

mash_input_data_loc <- snakemake@output[['mash_inputs']]

# subset to snps that were kept for all analyses
beta.hat <- read_tsv(beta_df_loc) %>%
  column_to_rownames("gv") %>% as.matrix

se.hat <- read_tsv(se_df_loc) %>%
  column_to_rownames("gv") %>% as.matrix

degf.hat <- read_tsv(sample_summary_loc) %>%
  filter(!dropped) %>%
  dplyr::count(type) %>%
  mutate(df=n-npcs-1) %>%
  select(c(type, n)) %>%
  column_to_rownames("type") %>% t %>%
  as_tibble() %>% 
  mutate(n=nrow(beta.hat)) %>% uncount(n) %>%
  mutate(gv=rownames(beta.hat)) %>%
  column_to_rownames("gv") %>% as.matrix

# band-aid - bug in mashr requires manual adjustment of standard errors
p2z = function(pval, Bhat) {
  z = abs(qnorm(pval / 2))
  z[which(Bhat < 0)] = -1 * z[which(Bhat < 0)]
  return(z)
}
## Shat = Bhat/Z where Z is the Z score corresponding to a p value from a t test done on (Bhat,Shat_orig,df)
z.adj = p2z(2 * pt(-abs(beta.hat/se.hat), degf.hat), beta.hat)
se.hat = beta.hat / z.adj

# Find the strongest test per gene
strongest_context_per_test <- tibble(gv=rownames(z.adj), 
                          max_value=apply(abs(z.adj), 1, max, na.rm=T),
                          max_context=colnames(z.adj)[max.col(replace(abs(z.adj), is.na(z.adj), -Inf))]) 

strongest_test_per_gene <- strongest_context_per_test %>%
  separate(gv, into=c("gene", "SNP"), sep="_", remove=F) %>%
  group_by(gene) %>% slice_max(max_value, n=1, with_ties=FALSE)

top.hits <- strongest_test_per_gene %>%
  pull(gv)

data.strong = mash_set_data(beta.hat[top.hits,], se.hat[top.hits,])

# Identify a random subset of 50K tests for model fitting, without too many missing values
keepers_cutoff <- intersect(rownames(beta.hat)[which(rowSums(!is.na(beta.hat)) >= n_contexts_cutoff)], 
                         rownames(se.hat)[which(rowSums(!is.na(se.hat)) >= n_contexts_cutoff)])

random.hits = sample(keepers_cutoff, min(50000, nrow(beta.hat)))

data.random = mash_set_data(beta.hat[random.hits,], se.hat[random.hits,])

# Get a stricter set of tests that were measured in every context (these will be used for learning regulatory patterns)
keepers_nonmissing <- intersect(rownames(beta.hat)[which(rowSums(is.na(beta.hat)) == 0)],
                                rownames(se.hat)[which(rowSums(is.na(se.hat)) == 0)]) 

strongest_tests_nonmissing <- strongest_context_per_test %>%
  filter(gv %in% keepers_nonmissing) %>%
  separate(gv, into=c("gene", "SNP"), sep="_", remove=F) %>%
  group_by(gene) %>% slice_max(max_value, n=1, with_ties=FALSE)
top.hits.nonmissing <- strongest_tests_nonmissing %>%
  arrange(desc(max_value)) %>%
  slice_head(n=2000) %>%
  pull(gv)

data.strong.nonmissing = mash_set_data(beta.hat[top.hits.nonmissing,], se.hat[top.hits.nonmissing,])

# Last but not least, get a sample of 25K tests (still no missingness) to estimate correlation structure
corr.subset = sample(keepers_nonmissing, 25000)

# estimate residual correlation
data.temp = mash_set_data(beta.hat[corr.subset,], se.hat[corr.subset,])
U.c = cov_canonical(data.temp, cov_methods = c("singletons", "equal_effects"))
V.em.full = mash_estimate_corr_em(data.temp, U.c, details = TRUE)
Vhat = V.em.full$V

# restore after df bug fix
# data.random = mash_set_data(beta.hat[random.hits,], se.hat[random.hits,], V=Vhat, df=degf.hat[random.hits,])
# data.strong = mash_set_data(beta.hat[top.hits,], se.hat[top.hits,], V=Vhat, df=degf.hat[top.hits,])

data.strong = mash_update_data(data.strong, V=Vhat)
data.random = mash_update_data(data.random, V=Vhat)
data.strong.nonmissing = mash_update_data(data.strong.nonmissing, V=Vhat)

# get data-driven covariances
# U.pca = cov_pca(data.strong.5k, 5)
# U.ed = cov_ed(data.strong, Ulist_init=U.pca)
# U.flash = cov_flash(data.strong, remove_singleton=T)
U.flash = cov_flash(data.strong.nonmissing, factors="nonneg", remove_singleton=T)

# fit mash model
U.c = cov_canonical(data.random, cov_methods = c("singletons", "equal_effects"))

save(data.temp, data.random, data.strong, data.strong.nonmissing, 
     strongest_test_per_gene, Vhat, V.em.full, U.flash, U.c, 
     file=mash_input_data_loc)
