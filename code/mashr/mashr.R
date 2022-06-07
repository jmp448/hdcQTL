library(tidyverse)
library(mashr)
library(vroom)
set.seed(1234)

qtl_files <- snakemake@input
trained_model_loc <- snakemake@output[['trained']]
tophit_fitted_loc <- snakemake@output[['tophits']]

print(qtl_files)

pull_type <- function(s) {
  str_split(s, "/")[[1]][5]
}

# compute standard errors, omitting tests where we can't obtain reasonable estimates
qtls <- vroom(qtl_files, id="path", col_select=c(SNP, gene, beta, `t-stat`, df, path)) %>%
  mutate(type=sapply(path, pull_type), .keep="unused") %>%
  mutate(se=if_else(as.logical((beta==0) * (`t-stat`==0)), rep(NA_real_, nrow(.)), beta/`t-stat`), .after='beta') %>%
  drop_na() %>%
  unite(gv, c(gene, SNP), sep="--", remove=FALSE)

print(table(qtl_files$type))

# 
# qtls <- map_dfr(qtl_files, load_qtls) %>%
#   mutate(se=if_else(`t-stat`==0, as.double(NA), beta/`t-stat`), .after='beta') %>%
#   mutate(se=if_else(se==0, as.double(NA), se)) %>%
#   drop_na() %>%
#   unite(gv, c(gene, SNP), sep="--", remove=FALSE)

# subset to snps that were kept for all analyses
keepers <- qtls %>%
  dplyr::count(gv) %>% 
  filter(n==length(qtl_files)) %>%
  pull(gv)

qtls <- qtls %>%
  filter(gv %in% keepers)

beta.hat <- qtls %>%
  select(gv, beta, type) %>%
  pivot_wider(names_from=type, values_from=beta, values_fill=NA) %>%
  column_to_rownames("gv") %>% as.matrix

se.hat <- qtls %>%
  select(gv, se, type) %>%
  pivot_wider(names_from=type, values_from=se, values_fill=NA) %>%
  column_to_rownames("gv") %>% as.matrix

df.hat <- qtls %>%
  select(gv, df, type) %>%
  pivot_wider(names_from=type, values_from=df, values_fill=NA) %>%
  column_to_rownames("gv") %>% as.matrix

# identify top hit per gene across all types/days
top.hits <- qtls %>%
  select(gv, gene, `t-stat`, type) %>%
  mutate(`t-stat`=abs(`t-stat`)) %>%
  group_by(gene) %>%
  slice_max(`t-stat`, n=1, with_ties=FALSE) %>%
  pull(gv)

# identify a random subset of 250K tests
random.hits = sample(qtls$gv, min(250000, length(unique(qtls$gv))))

# estimate null correlation
data.temp = mash_set_data(beta.hat[random.hits,],se.hat[random.hits,])
Vhat = estimate_null_correlation_simple(data.temp)

data.random = mash_set_data(beta.hat[random.hits,], se.hat[random.hits,], V=Vhat, df=df.hat[random.hits,])
data.strong = mash_set_data(beta.hat[top.hits,], se.hat[top.hits,], V=Vhat, df=df.hat[top.hits,])

# get data-driven covariances
U.pca = cov_pca(data.strong, 5)
U.ed = cov_ed(data.strong, U.pca)

# fit mash model
U.c = cov_canonical(data.random)
m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
saveRDS(m, trained_model_loc)

m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
saveRDS(m2, tophit_fitted_loc)