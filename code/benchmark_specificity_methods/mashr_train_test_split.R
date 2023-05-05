library(tidyverse)
library(mashr)
set.seed(1234)

# Load info from snakemake
beta_df_loc <- snakemake@input[['beta_df']]
se_df_loc <- snakemake@input[['se_df']]
gene_locs <- snakemake@input[['gene_locs']]

output_loc <- snakemake@output[[1]]

# Load data
beta.hat <- read_tsv(beta_df_loc) %>%
  column_to_rownames("gv") %>% as.matrix

se.hat <- read_tsv(se_df_loc) %>%
  column_to_rownames("gv") %>% as.matrix

# Filter to tests which have no NA values in either matrix
all_context_snps <- rownames(beta.hat)[which(rowSums(is.na(beta.hat)) + rowSums(is.na(se.hat)) == 0)]
beta.hat <- beta.hat[all_context_snps,]
se.hat <- se.hat[all_context_snps,]

# Divide up the previously used training data (200K variants) into training/ testing data sets
is.even <- function(chrom) {
  (as.numeric(str_replace(chrom, "chr", "")) %% 2) == 0
}

gene_assignments <- read_tsv(gene_locs) %>%
  mutate(even=is.even(chr)) %>%
  select(hgnc, even)

even.chroms <- tibble("gv"=rownames(beta.hat)) %>%
  tidyr::separate(gv, into=c("gene", "SNP"), sep="_") %>%
  left_join(gene_assignments, by=c("gene"="hgnc")) %>% 
  pull(even)

Bhat.train <- beta.hat[which(!even.chroms),]
Shat.train <- se.hat[which(!even.chroms),]
Bhat.test <- beta.hat[which(even.chroms),]
Shat.test <- se.hat[which(even.chroms),]

# Further reduce the train and test dataset sizes
train.subset <- sample(seq(1, nrow(Bhat.train)), 20000, replace=F)
test.subset <- sample(seq(1, nrow(Bhat.test)), 20000, replace=F)

Bhat.train <- Bhat.train[train.subset,]
Shat.train <- Shat.train[train.subset,]
Bhat.test <- Bhat.test[test.subset,]
Shat.test <- Shat.test[test.subset,]

save(Bhat.train, Shat.train, Bhat.test, Shat.test, file=output_loc)
