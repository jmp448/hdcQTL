### 
# Perform some EDA on pseudobulk samples to identify any outliers
# 
# Relies on having previously done aggregation in a notebook like `analysis/highpass/aggregate_*`
#
###

library(tidyverse)
library(Matrix)
library(Matrix.utils)
library(matrixStats)
set.seed(2021)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")
#use_condaenv("/project2/gilad/jpopp/ebQTL/.snakemake/conda/980130bc9d8d3f8069c208a84393fb39")

# Read in the input files from Snakefile
pseudobulk_loc <- snakemake@input[["pseudobulk"]]
sample_summary_loc <- as.character(snakemake@input[["sample_summary"]])
sample_summary_manual_loc <- as.character(snakemake@output[["sample_summary_manual"]])
table_prefix <- snakemake@params[['table_prefix']]
plot_prefix <- snakemake@params[['fig_prefix']]

# Get a list of cell types
cell.types <- read_tsv(pseudobulk_loc, n_max=0, col_select=-c(1)) %>%
  colnames() %>%
  sapply(function(s){str_split(s, "_")[[1]][2]}) %>% 
  unique()

# Normalize
invnorm_transform <- function(x) {
  qqnorm(x, plot.it=F)$x
}

# Read in pseudobulk data
pseudobulk <- read_tsv(pseudobulk_loc) %>% 
  arrange(gene)
type.ind <- colnames(pseudobulk)[-c(1)]

# Drop low-cell num samples
samples_subset <- read_tsv(sample_summary_loc) %>%
  filter(!dropped)

# Set up sample summary for manual check
if (!file.exists(sample_summary_manual_loc)) {
  file.copy(sample_summary_loc, sample_summary_manual_loc)
}

# Decompose expression into shared and cell type specific components
pseudobulk.decomp <- pseudobulk %>%
  column_to_rownames("gene") %>% t %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(cols=!sample, names_to="gene", values_to="expression") %>%
  mutate(ind=str_sub(sample, 1, 5)) %>%
  group_by(ind, gene) %>%
  mutate(shared=mean(expression)) %>%
  mutate(expression = expression - shared) %>%
  ungroup

pseudobulk.shared <- pseudobulk.decomp %>%
  select(ind, gene, shared) %>%
  distinct %>%
  mutate(ind=paste0(ind, "_Shared")) %>%
  rename(sample=ind, expression=shared)

pseudobulk.all <- pseudobulk.decomp %>%
  select(sample, gene, expression) %>%
  bind_rows(pseudobulk.shared) %>%
  arrange(sample) %>%
  pivot_wider(names_from=sample, values_from=expression)

# Prep for QTL calling
for (d in c(cell.types)) {
  dir.create(paste0(table_prefix, "/", d), recursive=T, showWarnings=FALSE)
  dir.create(paste0(plot_prefix, "/", d), recursive=T, showWarnings=FALSE)
  
  # get expression
  expression <- pseudobulk.all %>% 
    select(c(gene, ends_with(paste0("_", d)))) %>% 
    rename_with(function(x){if_else(x=="gene", true=x, false=paste0("NA", str_sub(x, 1, 5)))})
  expression.mat <- expression %>%
    column_to_rownames("gene") %>%
    as.matrix
  
  # gene filter - nonzero variance, nonzero expression in at least 5 samples
  gene_vars <- rowVars(expression.mat) > 0
  gene_nonzeros <- rowSums(expression.mat > 0) >= 5
  gene_keepers <- rownames(expression.mat)[gene_vars & gene_nonzeros]
  inds <- colnames(expression.mat)
  expression.mat <- expression.mat[gene_keepers,]
  
  # expression normalization
  expression.mat <- apply(expression.mat, 2, invnorm_transform)
  expression.mat <- t(apply(expression.mat, 1, invnorm_transform))
  rownames(expression.mat) <- gene_keepers
  colnames(expression.mat) <- colnames(expression)[-c(1)]
  expression <- as_tibble(expression.mat, rownames="gene")
  
  # expression PCs - produces a scree plot and a PC biplot, and saves the covariates
  pcomp <- regular.pca(expression)
  
  # viz PC biplot, scree, and num cells for the cell type
  png(file=paste0(plot_prefix, "/", d, "/qc.png"), width=480*3, height=480)
  par(mfrow=c(1,3))
  plot(pcomp$d ^ 2 / sum(pcomp$d ^ 2), main=d, xlab="PC", ylab="Variance Explained")
  plot(pcomp$u$PC1, pcomp$u$PC2, main=d, xlab="PC1", ylab="PC2")
  if (d != "Shared") {
    barplot(filter(samples_subset, type==d)$n_cells ~ as.factor(filter(samples_subset, type==d)$individual), 
            las=2, xlab="Individual", ylab="Number of Cells")
  }
  dev.off()
  
  covariates <- pcomp$u[,1:5] %>% 
    column_to_rownames("sample") %>% t %>% as_tibble(rownames="covariate") %>%
    write_tsv(paste0(table_prefix, "/", d, "/all_pcs.tsv"))
}