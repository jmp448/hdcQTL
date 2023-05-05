library(tidyverse)
library(mashr)
library(vroom)
set.seed(1234)

beta_df_loc <- snakemake@input[['beta_df']]
se_df_loc <- snakemake@input[['se_df']]
sample_summary_loc <- snakemake@input[['sample_summary']]

npcs <- as.numeric(snakemake@wildcards[['npcs']])

mash_input_data_loc <- snakemake@output[['mash_inputs']]

beta_df_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/tensorqtl_nominal.betas.tsv"
se_df_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/tensorqtl_nominal.standard_errors.tsv"
sample_summary_loc <- "data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/sample_summary_manual.tsv"

# subset to snps that were kept for all analyses
beta.hat <- read_tsv(beta_df_loc) %>%
  column_to_rownames("gv") %>% as.matrix

se.hat <- read_tsv(se_df_loc) %>%
  column_to_rownames("gv") %>% as.matrix

shared_snps <- intersect(rownames(beta.hat), rownames(se.hat))
beta.hat <- beta.hat[shared_snps,]
se.hat <- se.hat[shared_snps,]

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
se.hat[which(z.adj == 0)] = 10

# identify top hit per gene across all cell groups
top.hits <- as_tibble(z.adj, rownames="gv") %>%
  pivot_longer(!gv, names_to="type", values_to="z") %>%
  mutate(z=abs(z)) %>%
  separate(gv, into=c("gene", "SNP"), sep="_", remove=F) %>%
  group_by(gene) %>%
  slice_max(z, n=1, with_ties=FALSE) %>%
  pull(gv)

# identify a random subset of 250K tests for model fitting
random.hits = sample(rownames(beta.hat), min(250000, nrow(beta.hat)))

# identify an even smaller subset of 25K tests for estimating correlation
corr.subset = sample(random.hits, 25000)

# estimate residual correlation
data.temp = mash_set_data(beta.hat[corr.subset,], se.hat[corr.subset,])
U.c = cov_canonical(data.temp)
V.em.full = mash_estimate_corr_em(data.temp, U.c, details = TRUE)
Vhat = V.em.full$V

# restore after df bug fix
# data.random = mash_set_data(beta.hat[random.hits,], se.hat[random.hits,], V=Vhat, df=degf.hat[random.hits,])
# data.strong = mash_set_data(beta.hat[top.hits,], se.hat[top.hits,], V=Vhat, df=degf.hat[top.hits,])

data.random = mash_set_data(beta.hat[random.hits,], se.hat[random.hits,], V=Vhat)
data.strong = mash_set_data(beta.hat[top.hits,], se.hat[top.hits,], V=Vhat)

# get data-driven covariances
U.pca = cov_pca(data.strong, 5)
U.ed = cov_ed(data.strong, U.pca)

# fit mash model
U.c = cov_canonical(data.random, cov_methods = c("identity", "singletons", "equal_effects"))

save(data.temp, data.random, data.strong, Vhat, V.em.full, U.pca, U.ed, U.c, file=mash_input_data_loc)
