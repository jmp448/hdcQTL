library(tidyverse)
library(mashr)

# Load info from snakemake
# load(snakemake@input[['data']])
# gene_locs <- snakemake@input[['gene_locs']]
# output_loc <- snakemake@output

# Load info locally
load("results/static/highpass_cellid_all/pseudobulk-scran/basic/mashr.training_data.RData")
gene_locs <- "/project2/gilad/jpopp/ebQTL/data/gencode/gencode.hg38.filtered.tss.tsv"
output_loc <- "results/static/highpass_cellid_all/pseudobulk-scran/basic/mashr_alpha_eval.tsv"

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
U.c = cov_canonical(data.train)

logliks <- tibble("alpha"=c(0, 0.25, 0.5, 0.75, 1),
                  "loglik"=rep(0, 5))

# Try multiple possible values for alpha
for (i in seq(1, nrow(logliks))) {
  data.train.a = mash_set_data(Bhat.train, Shat.train, alpha=logliks$alpha[i])
  V.em.full = mash_estimate_corr_em(data.train.a, U.c, max_iter=10, details = TRUE)
  V.em = V.em.full$V
  m.Vem = V.em.full$mash.model
  
  data.test.a = mash_set_data(Bhat.test, Shat.test, alpha=logliks$alpha[i], V=V.em)
  m.Vem.test = mash(data.test.a, g=get_fitted_g(m.Vem), fixg=TRUE)
  logliks$loglik[i] = get_loglik(m.Vem.test)
  print(paste0("done with alpha = ", logliks$alpha[i]))
}

# Save outputs - train/test data, fitted mash models, 
write_tsv(logliks, output_loc)