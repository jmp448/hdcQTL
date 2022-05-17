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

# Read in the cell types and annotation information from Snakefile
extract_celltype <- function(s) {
  str_split(s, "/")[[1]] %>% rev %>% magrittr::extract2(2)
} 

cell.types <- snakemake@output %>%
  map_chr(extract_celltype) %>%
  unique

annotation_prefix <- sample_summary_loc %>%
  str_match(".*/") %>%
  as.character

#print(paste0("here is the snakemake output: \n", snakemake@output))
#print(paste0("\n \n here are the cell types: \n", cell.types))

# # PSEUDOBULK
invnorm_transform <- function(x) {
  qqnorm(x)$x
}

counts.type <- read_tsv(pseudobulk_loc) %>% 
  arrange(gene)
type.ind <- colnames(counts.type)[-c(1)]

samples_subset <- read_tsv(sample_summary_loc) %>%
  filter(!dropped)

# filter to genes with nonzero variance and expression in at least 5 samples, compute expression pcs, and save important information
# NA prepended to sample names in accordance with genotype data
for (d in cell.types) {
  dir.create(paste0(annotation_prefix, d), showWarnings = FALSE)
  expression <- counts.type %>% 
    select(c(gene, ends_with(paste0("_", d)))) %>% 
    rename_with(function(x){if_else(x=="gene", true=x, false=paste0("NA", str_sub(x, 1, 5)))})
  expression.mat <- expression %>%
    column_to_rownames("gene") %>%
    as.matrix
  gene_vars <- rowVars(expression.mat) > 0
  gene_nonzeros <- rowSums(expression.mat > 0) >= 5
  gene_keepers <- rownames(expression.mat)[gene_vars & gene_nonzeros]
  inds <- colnames(expression.mat)
  
  expression.mat <- expression.mat[gene_keepers,]
  expression.mat <- apply(expression.mat, 2, invnorm_transform)
  expression.mat <- t(apply(expression.mat, 1, invnorm_transform))
  rownames(expression.mat) <- gene_keepers
  colnames(expression.mat) <- colnames(expression)[-c(1)]
  
  expression <- as_tibble(expression.mat, rownames="gene") %>%
    write_tsv(paste0(annotation_prefix, d, "/expression.tsv"))
  
  pcomp <- regular.pca(expression)
  scree <- pcomp$d ^ 2 / sum(pcomp$d ^ 2) %>% 
    write(paste0(annotation_prefix, d, "/pve.tsv"), ncolumns=1)
  covariates <- pcomp$u[,1:5] %>% 
    column_to_rownames("sample") %>% t %>% as_tibble(rownames="covariate") %>%
    write_tsv(paste0(annotation_prefix, d, "/covariates.tsv"))
  individuals <- samples_subset %>%
    filter(type==!!d) %>%
    select(individual) %>%
    mutate(individual=paste0("NA", individual)) %>%
    write_tsv(paste0(annotation_prefix, d, "/individuals.tsv"), col_names=FALSE)
}