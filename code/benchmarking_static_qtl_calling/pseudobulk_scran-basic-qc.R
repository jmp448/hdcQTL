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

# Read in pseudobulk data & ensure it's sorted by gene and sample
pseudobulk <- read_tsv(pseudobulk_loc) %>% 
  arrange(gene)
sample_order <- pseudobulk %>%
  select(-c(gene)) %>%
  colnames %>%
  sort
pseudobulk <- relocate(pseudobulk, all_of(c("gene", sample_order)))

# Read in sample summary, sort and identify low-quality samples
samples_subset <- read_tsv(sample_summary_loc) %>%
  arrange(ind_type) %>%
  filter(!dropped)

# Set up sample summary for manual check
if (!file.exists(sample_summary_manual_loc)) {
  file.copy(sample_summary_loc, sample_summary_manual_loc)
}

# Prep for QTL calling
for (d in c(cell.types)) {
  dir.create(paste0(table_prefix, "/", d), recursive=T, showWarnings=FALSE)
  dir.create(paste0(plot_prefix, "/", d), recursive=T, showWarnings=FALSE)
  
  # get expression
  expression <- pseudobulk %>% 
    select(c(gene, ends_with(paste0("_", d)))) %>% 
    rename_with(function(x){if_else(x=="gene", true=x, false=paste0("NA", str_sub(x, 1, 5)))})
  expression.mat <- expression %>%
    column_to_rownames("gene") %>%
    as.matrix
  
  # gene filter - nonzero variance, nonzero median expression
  gene_vars <- rowVars(expression.mat) > 0
  gene_nonzeros <- rowMedians(expression.mat) > 0.25
  gene_keepers <- rownames(expression.mat)[gene_vars & gene_nonzeros]
  inds <- colnames(expression.mat)
  expression.mat <- expression.mat[gene_keepers,]
  
  # expression normalization
  # 1. normalize each gene to have zero mean, unit variance across individuals
  expression.mat <- apply(expression.mat, 1, invnorm_transform) %>% t
  # 2. re-insert gene names and convert to tibble
  rownames(expression.mat) <- gene_keepers
  expression <- as_tibble(expression.mat, rownames="gene")
  
  # expression PCs - produces a scree plot and a PC biplot, and saves the covariates
  pcomp <- regular.pca(expression)
  
  # viz PC biplot, scree, and num cells for the cell type
  png(file=paste0(plot_prefix, "/", d, "/qc.png"), width=480*3, height=480)
  par(mfrow=c(1,3))
  plot(pcomp$d ^ 2 / sum(pcomp$d ^ 2), main=d, xlab="PC", ylab="Variance Explained")
  plot(pcomp$u$PC1, pcomp$u$PC2, main=d, xlab="PC1", ylab="PC2")
  barplot(filter(samples_subset, type==d)$n_cells_filtered ~ as.factor(filter(samples_subset, type==d)$individual), 
          las=2, xlab="Individual", ylab="Number of Cells")
  dev.off()
  
  covariates <- pcomp$u %>% 
    column_to_rownames("sample") %>% t %>% as_tibble(rownames="covariate") %>%
    write_tsv(paste0(table_prefix, "/", d, "/all_pcs.tsv"))
}