library(tidyverse)
library(mashr)

# Load info from snakemake
# load(snakemake@input[['data']])
# gene_locs <- snakemake@input[['gene_locs']]
# output_loc <- snakemake@output

# Load info locally
load("results/static/highpass_cellid_all/pseudobulk-scran/basic/mashr.training_data.RData")
gene_locs <- "/project2/gilad/jpopp/ebQTL/data/gencode/gencode.hg38.filtered.tss.tsv"
output_loc <- "results/static/highpass_cellid_all/pseudobulk-scran/basic/mashr_corr_eval.RData"

# Divide up the previously used training data (200K variants) into training/ testing data sets
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

# Further reduce the train and test dataset sizes
train.subset <- sample(seq(1, nrow(Bhat.train)), 20000)
test.subset <- sample(seq(1, nrow(Bhat.test)), 20000)
Bhat.train <- Bhat.train[train.subset,]
Shat.train <- Shat.train[train.subset,]
Bhat.test <- Bhat.test[test.subset,]
Shat.test <- Shat.test[test.subset,]

data.train <- mash_set_data(Bhat.train, Shat.train)
data.test <- mash_set_data(Bhat.test, Shat.test)

# Estimate correlation using the basic approach
V.simple = estimate_null_correlation_simple(data.train)
data.train.Vsimple = mash_update_data(data.train, V=V.simple)
data.test.Vsimple = mash_update_data(data.test, V=V.simple)

U.c = cov_canonical(data.train)

# Estimate mash parameters using basic approach (with and without simple correlation estimate)
m.orig = mash(data.train, U.c)
m.Vsimple = mash(data.train.Vsimple, U.c)

# Estimate correlation & mash parameters using EM approach
V.em.full = mash_estimate_corr_em(data.train, U.c, details = TRUE)
V.em = V.em.full$V
m.Vem = V.em.full$mash.model
data.test.Vem = mash_update_data(data.test, V=V.em)

# Compare model fits on held-out test data
m.orig.test = mash(data.test, g=get_fitted_g(m.orig), fixg=TRUE)
m.Vsimple.test = mash(data.test.Vsimple, g=get_fitted_g(m.Vsimple), fixg=TRUE)
m.Vem.test = mash(data.test.Vem, g=get_fitted_g(m.Vem), fixg=TRUE)

# Save outputs - train/test data, fitted mash models, 
save(Bhat.train, Shat.train,
     Bhat.test, Shat.test,
     V.simple, V.em,
     m.orig.test, m.Vsimple.test, m.Vem.test,
     file=output_loc)