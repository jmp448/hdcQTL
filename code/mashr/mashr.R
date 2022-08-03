library(tidyverse)
library(mashr)
library(vroom)
set.seed(1234)

# cell.types <- read_tsv("data/single_cell_objects/highpass_cellid_all.pseudobulk-scran.tsv", n_max=0, col_select=-c(gene)) %>%
#   colnames %>%
#   sapply(function(s){str_sub(s, start=7)}) %>%
#   unique()
# qtl_files <- paste0("results/static/highpass_cellid_all/pseudobulk-scran/basic/", cell.types, "/matrixeqtl.cis_qtl_pairs.all.mtc.tsv")

qtl_files <- snakemake@input
trained_model_loc <- snakemake@output[['trained']]
tophit_fitted_loc <- snakemake@output[['tophits']]
random_hits_loc <- snakemake@output[['random']]
top_hits_loc <- snakemake@output[['top']]
mashr_input_data_loc <- snakemake@output[['datasets']]
qtl_master_loc <- snakemake@output[['all_qtls']]

pull_type <- function(s) {
  str_split(s, "/")[[1]][6]
}

# compute standard errors, omitting tests where we can't obtain reasonable estimates
qtls <- vroom(qtl_files, id="path", col_select=c(SNP, gene, beta, `t-stat`, df, path)) %>%
  mutate(type=sapply(path, pull_type), .keep="unused") %>%
  mutate(se=if_else(as.logical((beta==0) * (`t-stat`==0)), rep(NA_real_, nrow(.)), beta/`t-stat`), .after='beta') %>%
  drop_na() %>%
  unite(gv, c(gene, SNP), sep="--", remove=FALSE)

saveRDS(qtls, qtl_master_loc)

# we'll use this later, mapping each cell type to the degrees of freedom
df.type <- qtls %>%
  select(c(type, df)) %>%
  distinct()

# subset to snps that were kept for all analyses
keepers <- qtls %>%
  dplyr::count(gv) %>% 
  filter(n >= 20) %>%
  pull(gv)

qtls <- qtls %>%
  filter(gv %in% keepers)

beta.hat <- qtls %>%
  select(gv, beta, type) %>%
  pivot_wider(names_from=type, values_from=beta, values_fill=0) %>%
  column_to_rownames("gv") %>% as.matrix

se.hat <- qtls %>%
  select(gv, se, type) %>%
  pivot_wider(names_from=type, values_from=se, values_fill=10) %>%
  column_to_rownames("gv") %>% as.matrix

df.hat <- df.type %>%
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
z.adj = p2z(2 * pt(-abs(beta.hat/se.hat), df.hat), beta.hat)
se.hat = beta.hat / z.adj
se.hat[which(z.adj == 0)] = 10

# identify top hit per gene across all types/days
top.hits <- qtls %>%
  select(gv, gene, `t-stat`, type) %>%
  mutate(`t-stat`=abs(`t-stat`)) %>%
  group_by(gene) %>%
  slice_max(`t-stat`, n=1, with_ties=FALSE) %>%
  pull(gv)
write_tsv(tibble(top.hits), top_hits_loc, col_names=F)

# identify a random subset of 250K tests
random.hits = sample(qtls$gv, min(250000, length(unique(qtls$gv))))
write_tsv(tibble(random.hits), random_hits_loc, col_names=F)

# estimate null correlation
data.temp = mash_set_data(beta.hat[random.hits,],
                          se.hat[random.hits,])
Vhat = estimate_null_correlation_simple(data.temp)

# restore after df bug fix
# data.random = mash_set_data(beta.hat[random.hits,], se.hat[random.hits,], V=Vhat, df=df.hat[random.hits,])
# data.strong = mash_set_data(beta.hat[top.hits,], se.hat[top.hits,], V=Vhat, df=df.hat[top.hits,])

data.random = mash_set_data(beta.hat[random.hits,], se.hat[random.hits,], V=Vhat)
data.strong = mash_set_data(beta.hat[top.hits,], se.hat[top.hits,], V=Vhat)

save(data.temp, Vhat, data.random, data.strong, file=mashr_input_data_loc)

# get data-driven covariances
U.pca = cov_pca(data.strong, 5)
U.ed = cov_ed(data.strong, U.pca)

# fit mash model
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
saveRDS(m, trained_model_loc)

m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
saveRDS(m2, tophit_fitted_loc)