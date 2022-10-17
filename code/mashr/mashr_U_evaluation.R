library(tidyverse)
library(mashr)

# Load info from snakemake
# load(snakemake@input[['data']])
# gene_locs <- snakemake@input[['gene_locs']]
# output_loc <- snakemake@output

# Load info locally
load("results/static/highpass_cellid_all/pseudobulk-scran/basic/mashr.training_data.RData")
gene_locs <- "/project2/gilad/jpopp/ebQTL/data/gencode/gencode.hg38.filtered.tss.tsv"
output_loc <- "results/static/highpass_cellid_all/pseudobulk-scran/basic/mashr_U_eval.RData"

## Divide up the previously used training data (200K variants) into training/ testing data sets
is.even <- function(chrom) {
  (as.numeric(str_replace(chrom, "chr", "")) %% 2) == 0
}
gene_assignments <- read_tsv(gene_locs) %>%
  mutate(even=is.even(chr)) %>%
  select(hgnc, even)

even.chroms <- tibble(gv=rownames(data.random$Bhat)) %>%
  separate(gv, into=c("gene", "SNP"), sep="--") %>%
  left_join(gene_assignments, by=c("gene"="hgnc")) %>% 
  pull(even)

Bhat.train <- data.random$Bhat[!even.chroms,]
Shat.train <- data.random$Shat[!even.chroms,]
Bhat.test <- data.random$Bhat[even.chroms,]
Shat.test <- data.random$Shat[even.chroms,]

## Further reduce the train and test dataset sizes
train.subset <- sample(seq(1, nrow(Bhat.train)), 20000)
test.subset <- sample(seq(1, nrow(Bhat.test)), 20000)
Bhat.train <- Bhat.train[train.subset,]
Shat.train <- Shat.train[train.subset,]
Bhat.test <- Bhat.test[test.subset,]
Shat.test <- Shat.test[test.subset,]

data.train <- mash_set_data(Bhat.train, Shat.train)
data.test <- mash_set_data(Bhat.test, Shat.test)

## Get the strong tests from just the odd chromosomes for learning data-driven covariance / training
even.chroms.strong <- tibble(gv=rownames(data.strong$Bhat)) %>%
  separate(gv, into=c("gene", "SNP"), sep="--") %>%
  left_join(gene_assignments, by=c("gene"="hgnc")) %>% 
  pull(even)

data.train.strong <- mash_set_data(data.strong$Bhat[!even.chroms.strong,],
                                   data.strong$Shat[!even.chroms.strong,])

# Estimate correlation using the basic approach
V.simple = estimate_null_correlation_simple(data.train)
data.train = mash_update_data(data.train, V=V.simple)
data.train.strong = mash_update_data(data.train.strong, V=V.simple)
data.test = mash_update_data(data.test, V=V.simple)

U.c = cov_canonical(data.train, cov_methods = c("identity", "singletons", "equal_effects"))

# Estimate mash parameters without data-driven covariance matrices
m.simple = mash(data.train, U.c, outputlevel = 1)

# Estimate mash parameters with default data-driven covariance matrices
U.pca.default = cov_pca(data.train.strong, 5)
U.ed.default = cov_ed(data.train.strong, U.pca.default)
m.default = mash(data.random, Ulist = c(U.ed.default, U.c), outputlevel = 1)

# Estimate mash parameters with data-driven covariance matrices similar to the GTEx paper
U.pca.paper = cov_pca(data.train.strong, 3)[4]
U.flash = cov_flash(data.train.strong, Kmax=5)
U_x = apply(data.train.strong$Bhat, 2, function(x) x - mean(x))
U_xx = t(U_x) %*% U_x / nrow(U_x)

U.ed.paper = cov_ed(data.train.strong, Ulist_init=c(U.pca.paper, U.flash[1], list(xx=U_xx)))
U.sfa.paper = U.flash[2:6]
m.paper = mash(data.random, Ulist = c(U.ed.paper, U.sfa.paper, U.c), outputlevel = 1)

# Estimate mash parameters with default data-driven covariance matrices but for z-stats
data.train.strong.EZ = mash_set_data(data.strong$Bhat[!even.chroms.strong,],
                                     data.strong$Shat[!even.chroms.strong,],
                                     alpha=1, V=V.simple)
U.pca.default.EZ = cov_pca(data.train.strong.EZ, 5)
U.ed.default.EZ = cov_ed(data.train.strong.EZ, U.pca.default.EZ)
m.default.EZ = mash(data.random, Ulist = c(U.ed.default.EZ, U.c), outputlevel = 1)

# Estimate mash parameters with data-driven covariance matrices similar to the GTEx paper but again with z-stats
U.pca.paper.EZ = cov_pca(data.train.strong.EZ, 3)[4]
U.flash.EZ = cov_flash(data.train.strong.EZ, Kmax=5)
U_x.EZ = apply(data.train.strong.EZ$Bhat, 2, function(x) x - mean(x))
U_xx.EZ = t(U_x.EZ) %*% U_x.EZ / nrow(U_x.EZ)

U.ed.paper.EZ = cov_ed(data.train.strong.EZ, Ulist_init=c(U.pca.paper.EZ, U.flash.EZ[1], list(xx=U_xx.EZ)))
U.sfa.paper.EZ = U.flash.EZ[2:6]
m.paper.EZ = mash(data.random, Ulist = c(U.ed.paper.EZ, U.sfa.paper.EZ, U.c), outputlevel = 1)

# Compare model fits on held-out test data
m.simple.test = mash(data.test, g=get_fitted_g(m.simple), fixg=TRUE)
m.default.test = mash(data.test, g=get_fitted_g(m.default), fixg=TRUE)
m.paper.test = mash(data.test, g=get_fitted_g(m.paper), fixg=TRUE)
m.default.test.EZ = mash(data.test, g=get_fitted_g(m.default.EZ), fixg=TRUE)
m.paper.test.EZ = mash(data.test, g=get_fitted_g(m.paper.EZ), fixg=TRUE)

# Save outputs - train/test data, fitted mash models, 
save(m.simple.test,
     m.default.test, m.paper.test,
     m.default.test.EZ, m.paper.test.EZ,
     file=output_loc)