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
library(edgeR)
set.seed(2021)
#use_condaenv("/project2/gilad/jpopp/ebQTL/.snakemake/conda/980130bc9d8d3f8069c208a84393fb39")

# Read in the input files from Snakefile
pseudobulk_loc <- snakemake@input[["pseudobulk"]]
sample_summary_loc <- as.character(snakemake@input[["sample_summary_manual"]])
celltypes_loc <- as.character(snakemake@input[["celltypes"]])
table_prefix <- snakemake@params[['table_prefix']]
plot_prefix <- snakemake@params[['fig_prefix']]
all_expression_loc <- snakemake@output[["all_expression"]]

# Get a list of cell types
cell.types <- read_tsv(celltypes_loc, n_max=0, col_select=-c(1)) %>%
  colnames() %>%
  sapply(function(s){str_split(s, "_")[[1]][2]}) %>% 
  unique()

# Normalize
invnorm_transform <- function(x) {
  qqnorm(x, plot.it=F)$x
}

scale_tidy <- function(x, scale=TRUE) {
  # a version of the scale function that doesn't expect a matrix
  if (scale==TRUE) {
    (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
  } else {
    (x - mean(x, na.rm=TRUE))
  }
}

# Read in pseudobulk data
pseudobulk <- read_tsv(pseudobulk_loc) %>% 
  arrange(gene)
type.ind <- colnames(pseudobulk)[-c(1)]

# Identify low-quality samples
samples_subset <- read_tsv(sample_summary_loc) %>%
  filter(!dropped)

# In each context, filter to genes with nonzero variance and median expression 
# over ten, then apply TMM normalization, then re-merge all cell types together
pseudobulk_tmm <- tibble(gene=pseudobulk$gene)
for (d in c(cell.types)) {
  # get expression
  expression <- pseudobulk %>% 
    select(c(gene, ends_with(paste0("_", d)))) %>% 
    rename_with(function(x){if_else(x=="gene", true=x, false=paste0("NA", str_sub(x, 1, 5)))})
  expression.mat <- expression %>%
    column_to_rownames("gene") %>%
    as.matrix
  
  # gene filter - nonzero variance, nonzero median expression
  gene_vars <- rowVars(expression.mat) > 0
  gene_nonzeros <- rowMedians(expression.mat) > 10
  gene_keepers <- rownames(expression.mat)[gene_vars & gene_nonzeros]
  inds <- colnames(expression.mat)
  expression.mat <- expression.mat[gene_keepers,]
  
  # TMM normalization to get log TMM-normalized counts
  expression.mat <- DGEList(counts=expression.mat) %>% calcNormFactors(method="TMM") %>% cpm(log=T)
  
  expression <- as_tibble(expression.mat, rownames="gene") %>%
    rename_with(function(x){if_else(x=="gene", true=x, false=paste0(str_sub(x, -5), "_", d))})
  
  #pseudobulk_tmm <- left_join(pseudobulk_tmm, expression, by="gene")
  pseudobulk_tmm <- inner_join(pseudobulk_tmm, expression, by="gene")
}

# Make this a long matrix to aid in decomposition
pseudobulk_tmm <- pseudobulk_tmm %>%
  pivot_longer(cols=!gene, names_to="sample", values_to="expression") %>%
  mutate(ind=str_sub(sample, 1, 5), type=str_sub(sample, 7))

if (snakemake@wildcards[['decomp']] == "fastgxc") {
  pseudobulk <- ungroup(pseudobulk_tmm)
} else if (snakemake@wildcards[['decomp']] == "fastgxc_context_centered") {
  pseudobulk <- pseudobulk_tmm %>%
    mutate(expression=scale_tidy(expression, scale=FALSE)) %>%
    ungroup
} else if (snakemake@wildcards[['decomp']] == "fastgxc_context_standardized") {
  pseudobulk <- pseudobulk_tmm %>%
    mutate(expression=scale_tidy(expression, scale=TRUE)) %>%
    ungroup
} else {
  stop("invalid decomp wildcard")
}

# Decompose expression into shared and cell type specific components
pseudobulk.decomp <- pseudobulk %>%
  group_by(ind, gene) %>%
  mutate(shared=mean(expression)) %>%
  mutate(expression = expression - shared) %>%
  ungroup

pseudobulk.shared <- pseudobulk.decomp %>%
  select(ind, gene, shared) %>%
  distinct %>%
  mutate(sample=paste0(ind, "_Shared"), .keep="unused") %>%
  rename(expression=shared)

pseudobulk.all <- pseudobulk.decomp %>%
  select(sample, gene, expression) %>%
  bind_rows(pseudobulk.shared) %>%
  arrange(sample) %>%
  pivot_wider(id_cols=gene, names_from="sample", values_from="expression") %>%
  write_tsv(all_expression_loc)

for (d in c(cell.types, "Shared")) {
  dir.create(paste0(table_prefix, "/", d), recursive=T, showWarnings=FALSE)
  dir.create(paste0(plot_prefix, "/", d), recursive=T, showWarnings=FALSE)
  
  expression <- pseudobulk.all %>% 
    select(c(gene, ends_with(paste0("_", d)))) %>%
    rename_with(function(x){if_else(x=="gene", true=x, false=paste0("NA", str_sub(x, 1, 5)))})
  expression.mat <- expression %>%
    column_to_rownames("gene") %>%
    as.matrix
  
  # expression normalization
  gene_names <- rownames(expression.mat) 
  # normalize each gene to have zero mean, unit variance across individuals
  expression.mat <- t(scale(t(expression.mat), center=T, scale=T))
  # quantile normalize within each individual across all genes to standard normal
  expression.mat <- apply(expression.mat, 2, invnorm_transform)
  # re-insert gene names and convert to tibble
  rownames(expression.mat) <- gene_names
  expression <- as_tibble(expression.mat, rownames="gene") %>%
    write_tsv(paste0(table_prefix, "/", d, "/expression.tsv"))
  
  # expression PCs - produces a scree plot and a PC biplot, and saves the covariates
  pcomp <- princomp(expression.mat)
  
  # viz PC biplot, scree, and num cells for the cell type
  png(file=paste0(plot_prefix, "/", d, "/qc_manual.png"), width=480*3, height=480)
  par(mfrow=c(1,3))
  plot(pcomp$sdev ^ 2 / sum(pcomp$sdev ^ 2), main=d, xlab="PC", ylab="Variance Explained")
  plot(pcomp$loadings[,"Comp.1"], pcomp$loadings[,"Comp.2"], main=d, xlab="PC1", ylab="PC2")
  if (d != "Shared") {
    barplot(filter(samples_subset, type==d)$n_cells_filtered ~ as.factor(filter(samples_subset, type==d)$individual), 
            las=2, xlab="Individual", ylab="Number of Cells")
  }
  dev.off()
  
  covariates <- pcomp$loadings[,1:5] %>% t %>% 
    as_tibble(rownames="covariate") %>%
    mutate(covariate=str_replace(covariate, "Comp.", "PC")) %>%
    write_tsv(paste0(table_prefix, "/", d, "/covariates.tsv"))
  
  # NA prepended to sample names in accordance with genotype data
  if (d == "Shared") {
    individuals <- tibble(individual=unique(samples_subset$individual)) %>%
      mutate(individual=paste0("NA", individual)) %>%
      write_tsv(paste0(table_prefix, "/Shared/individuals.tsv"), col_names=FALSE)
  } else {
    individuals <- samples_subset %>%
      filter(type==!!d) %>%
      mutate(individual=paste0("NA", individual), .keep="used") %>%
      write_tsv(paste0(table_prefix, "/", d, "/individuals.tsv"), col_names=FALSE)
  }
}