### 
# Prep for Matrix eQTL by measuring expression covariates, etc.
# 
# Annotation prefix allows for multiple approaches to cell type annotation to be efficiently used in the workflow (ie multiple resolution clusterings)
# Relies on having previously done aggregation in a notebook like `analysis/highpass/aggregation_*`
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
celltype_summary_loc <- as.character(snakemake@input[["celltype_summary"]])
table_prefix <- snakemake@params[['table_prefix']]
plot_prefix <- snakemake@params[['fig_prefix']]

# Read in the output files from Snakefile
extract_celltype <- function(s) {
  str_split(s, "/")[[1]] %>% rev %>% magrittr::extract2(2)
} 
cell.types <- snakemake@output[["tables"]] %>%
  map_chr(extract_celltype) %>%
  unique

print(paste0("here is the snakemake output: ", snakemake@output[["tables"]]))
print(paste0("here are the cell types: ", cell.types))

# PSEUDOBULK
invnorm_transform <- function(x) {
  qqnorm(x, plot.it=F)$x
}

counts.type <- read_tsv(pseudobulk_loc) %>% 
  arrange(gene)
type.ind <- colnames(counts.type)[-c(1)]

samples_subset <- read_tsv(sample_summary_loc) %>%
  filter(!dropped)

for (d in cell.types) {
  dir.create(paste0(table_prefix, "/", d), recursive=T, showWarnings=FALSE)
  dir.create(paste0(plot_prefix, "/", d), recursive=T, showWarnings=FALSE)
  
  expression <- counts.type %>% 
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
  expression <- as_tibble(expression.mat, rownames="gene") %>%
    write_tsv(paste0(table_prefix, "/", d, "/expression.tsv"))
  
  # expression PCs - produces a scree plot and a PC biplot, and saves the covariates
  pcomp <- regular.pca(expression)
  
  png(file=paste0(plot_prefix, "/", d, "/scree.png"))
  plot(pcomp$d ^ 2 / sum(pcomp$d ^ 2), main=d, xlab="PC", ylab="Variance Explained")
  dev.off()
  
  png(file=paste0(plot_prefix, "/", d, "/pcs.png"))
  plot(pcomp$u$PC1, pcomp$u$PC2, main=d, xlab="PC1", ylab="PC2")
  dev.off()
  
  covariates <- pcomp$u[,1:5] %>% 
    column_to_rownames("sample") %>% t %>% as_tibble(rownames="covariate") %>%
    write_tsv(paste0(table_prefix, "/", d, "/covariates.tsv"))
  
  # NA prepended to sample names in accordance with genotype data
  individuals <- samples_subset %>%
    filter(type==!!d) %>%
    mutate(individual=paste0("NA", individual), .keep="used") %>%
    write_tsv(paste0(table_prefix, "/", d, "/individuals.tsv"), col_names=FALSE)
}