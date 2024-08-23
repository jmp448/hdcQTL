### 
# Get pseudobulk aggregated across all donors, for visualizing expression trends more easily
###

library(tidyverse)
library(Matrix)
library(matrixStats)
library(edgeR)
set.seed(2021)
source("code/dynamic_qtl_calling/clpca.R")

# raw_expression_loc <- "data/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/pseudobulk_raw.tsv"
# norm_expression_loc <- "data/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/nipals/pseudobulk_normalized_alldonors.tsv"
# sample_metadata_loc <- "/project2/gilad/katie/ebQTL/CombinedFormationAndCollectionMetadata_102andPilot_SWAPSANDCONTAMINATIONADDED_012522.csv"
# pseudobulk_loc <- "data/dynamic_qtl_calling/eb-neur_15binstrimmed/pseudobulk_tmm/eb-neur_15binstrimmed.pseudobulk_tmm.tsv"
# sample_summary_loc <- "data/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/sample_summary.tsv"

raw_expression_loc <- snakemake@output[["raw_expression"]]
norm_expression_loc <- snakemake@output[["norm_expression"]]

# Inverse normal function
invnorm_transform <- function(x) {
  qqnorm(x, plot.it=F)$x
}

expression_raw <- read_tsv(raw_expression_loc) %>% 
  gather(key="donor_bin", value="counts", -gene) %>%
  separate(donor_bin, into=c("donor", "bin")) %>% 
  group_by(gene, bin) %>%
  summarize(sum_counts=sum(counts)) %>%
  spread(key=bin, value=sum_counts)

expression.mat <- column_to_rownames(expression_raw, "gene")

# expression normalization
gene_names <- expression_raw$gene
bin_indices <- colnames(expression_raw)[-c(1)]
# 1. TMM normalization to get log TMM-normalized counts
expression.mat <- DGEList(counts=expression.mat) %>% calcNormFactors(method="TMM") %>% cpm(log=T)
# 2. normalize each gene to have zero mean, unit variance across individuals
expression.mat <- apply(expression.mat, 1, invnorm_transform) %>% t
# 3. re-insert gene names and convert to tibble
rownames(expression.mat) <- gene_names
colnames(expression.mat) <- bin_indices

expression_norm <- as_tibble(expression.mat, rownames="gene") %>%
  write_tsv(norm_expression_loc)
