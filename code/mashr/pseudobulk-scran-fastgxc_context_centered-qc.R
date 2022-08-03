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

# Identify low-quality samples
samples_subset <- read_tsv(sample_summary_loc) %>%
  filter(!dropped)

# Set up sample summary for manual check
if (!file.exists(sample_summary_manual_loc)) {
  file.copy(sample_summary_loc, sample_summary_manual_loc)
}

# First, filter in each context to genes with nonzero variance 
#   and at least 5 samples with nonzero mean expression
pseudobulk <- pseudobulk %>%
  pivot_longer(cols=!gene, names_to="sample", values_to="expression") %>%
  filter(sample %in% samples_subset$ind_type) %>%
  mutate(ind=str_sub(sample, 1, 5), type=str_sub(sample, 7)) %>%
  group_by(type, gene) %>% 
  mutate(nonzero_exp=sum(expression > 0), nonzero_var=var(expression)>0) %>%
  filter((nonzero_exp > 5) & (nonzero_var == T)) %>%
  select(-c(nonzero_exp, nonzero_var, type))

if (snakemake@wildcards[['decomp']] == "fastgxc") {
  pseudobulk <- ungroup(pseudobulk)
} else if (snakemake@wildcards[['decomp']] == "fastgxc_context_centered") {
  pseudobulk <- pseudobulk %>%
    mutate(expression=scale(expression, scale=FALSE)) %>%
    ungroup
} else if (snakemake@wildcards[['decomp']] == "fastgxc_context_standardized") {
  pseudobulk <- pseudobulk %>%
    mutate(expression=scale(expression, scale=FALSE)) %>%
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
  arrange(sample)

# Prep for QTL calling
for (d in c(cell.types, "Shared")) {
  dir.create(paste0(table_prefix, "/", d), recursive=T, showWarnings=FALSE)
  dir.create(paste0(plot_prefix, "/", d), recursive=T, showWarnings=FALSE)
  print(paste0("Reached ", d))
  
  # get expression
  expression <- pseudobulk.all %>% 
    filter(str_extract(sample, "[^_]+$")== !!d) %>%
    pivot_wider(names_from=sample, values_from=expression) %>%
    rename_with(function(x){if_else(x=="gene", true=x, false=paste0("NA", str_sub(x, 1, 5)))})
  if (d == "Shared") {
    expression <- drop_na(expression)
  }
  expression.mat <- expression %>%
    column_to_rownames("gene") %>%
    as.matrix
  
  # expression normalization
  gene_names <- rownames(expression.mat)
  expression.mat <- apply(expression.mat, 2, invnorm_transform)
  expression.mat <- t(apply(expression.mat, 1, invnorm_transform))
  rownames(expression.mat) <- gene_names
  colnames(expression.mat) <- colnames(expression)[-c(1)]
  expression <- as_tibble(expression.mat, rownames="gene")
  
  # expression PCs - produces a scree plot and a PC biplot, and saves the covariates
  pcomp <- princomp(expression.mat)
  
  # viz PC biplot, scree, and num cells for the cell type
  png(file=paste0(plot_prefix, "/", d, "/qc.png"), width=480*3, height=480)
  par(mfrow=c(1,3))
  plot(pcomp$sdev ^ 2 / sum(pcomp$sdev ^ 2), main=d, xlab="PC", ylab="Variance Explained")
  plot(pcomp$loadings[,"Comp.1"], pcomp$loadings[,"Comp.2"], main=d, xlab="PC1", ylab="PC2")
  if (d != "Shared") {
    barplot(filter(samples_subset, type==d)$n_cells ~ as.factor(filter(samples_subset, type==d)$individual), 
            las=2, xlab="Individual", ylab="Number of Cells")
  }
  dev.off()
  
  covariates <- pcomp$loadings[,1:5] %>% t %>% 
    as_tibble(rownames="covariate") %>%
    mutate(covariate=str_replace(covariate, "Comp.", "PC")) %>%
    write_tsv(paste0(table_prefix, "/", d, "/all_pcs.tsv"))
}