library(sva)
library(readr)
library(stringr)
library(dplyr)
library(vroom)
library(tidyr)
library(tibble)
library(purrr)
library(matrixStats)

set.seed(1234)

genotypes_loc <- snakemake@input[['genotypes']]
covariates_loc <- snakemake@input[['cov']]
expression_loc <- snakemake@input[['exp']]
sv_loc <- snakemake@output[['svs']]

wrangle_mod <- function(v, geno, cov, intercept=F) {
  g = dplyr::filter(geno, snp==!!v) %>% 
    rename(ID=snp)
  mod = bind_rows(cov, g) %>%
    column_to_rownames("ID") %>% t 
  mod
}

wrangle_null_mod <- function(cov, intercept=F) {
  null_mod = cov %>%
    column_to_rownames("ID") %>% t 
  null_mod
}

estimate_n_svs <- function(v, geno, cov, exp) {
  mod = wrangle_mod(v, geno, cov)
  num.sv(exp, mod)
}

compute_svs <- function(v, geno, cov, exp, n.sv) {
  mod = wrangle_mod(v, geno, cov)
  mod0 = wrangle_null_mod(cov)
  svs = sva(expression, mod, mod0, n.sv, method="irw")
  snp_sv = t(svs$sv)
  colnames(snp_sv) = colnames(geno)[2:ncol(geno)]
  sv_out = as_tibble(snp_sv, rownames="ID") %>%
    mutate(ID=paste0("SV", ID))
  geno_pc_out = mutate(cov, ID=paste0("GENO_", ID))
  cov_out = bind_rows(sv_out, geno_pc_out) %>%
    mutate(snp = v, .before=1)
  cov_out
}

# Load data
### Genotypes
print(paste0("genotypes at ", genotypes_loc))
genotypes <- vroom(genotypes_loc)

### Genotype PCs (known covariates)
print(paste0("covariates at ", covariates_loc))
covariates <- vroom(covariates_loc) %>%
  dplyr::filter(str_sub(ID, 1, 2) == "PC")

### Expression
print(paste0("expression at ", expression_loc))
expression <- vroom(expression_loc) %>%
  column_to_rownames("gene_id") %>%
  as.matrix
  
# Get a sample of SNPs to figure out the appropriate number of SVs to use in this tissue
sampled_snps <- sample(genotypes$snp, 10)
n_svs <- sapply(sampled_snps, estimate_n_svs, geno=genotypes, cov=covariates, exp=expression)
n_svs_tissue <- median(n_svs)

# Compute supervised SVs for each SNP
tissue_svs <- bind_rows(lapply(genotypes$snp, compute_svs, geno=genotypes, cov=covariates, exp=expression, n.sv=n_svs_tissue))
write_tsv(tissue_svs, sv_loc)