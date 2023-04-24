### 
# Aggregate pseudobulk within each cell type
# 
# Annotation prefix allows for multiple approaches to cell type annotation to be efficiently used in the workflow (ie multiple resolution clusterings)
# Uses TMM-normalization to account for (individual, cell type) samples having different depths
# Removes any samples with less than 10 (or `min_cells_per_sample` value) cells
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

# Identify low-quality samples & sort sample summary
samples_subset <- read_tsv(sample_summary_loc) %>%
  arrange(ind_type) %>%
  filter(!dropped)

# Set up sample summary for manual check
if (!file.exists(sample_summary_manual_loc)) {
  file.copy(sample_summary_loc, sample_summary_manual_loc)
}

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

### How many cell types is each gene measured in?
# pseudobulk_tmm %>%
#   drop_na() %>%
#   select(gene, type) %>%
#   distinct() %>%
#   count(gene) %>%
#   ggplot(aes(x=n)) +
#   geom_histogram(bins=5)

if (snakemake@wildcards[['decomp']] == "fastgxc") {
  pseudobulk <- ungroup(pseudobulk_tmm)
} else if (snakemake@wildcards[['decomp']] == "fastgxc_context_centered") {
  pseudobulk <- pseudobulk_tmm %>%
    group_by(gene, type) %>%
    mutate(expression=scale_tidy(expression, scale=FALSE)) %>%
    ungroup
} else if (snakemake@wildcards[['decomp']] == "fastgxc_context_standardized") {
  pseudobulk <- pseudobulk_tmm %>%
    group_by(gene, type) %>%
    mutate(expression=scale_tidy(expression, scale=TRUE)) %>%
    ungroup
} else {
  stop("invalid decomp wildcard")
}

# Decompose expression into shared and cell type specific components
pseudobulk.decomp <- pseudobulk %>%
  group_by(ind, gene) %>%
  mutate(shared=mean(expression, na.rm=T)) %>%
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
  
  # get expression
  expression <- pseudobulk.all %>% 
    filter(str_extract(sample, "[^_]+$")== !!d) %>%
    drop_na() %>%
    pivot_wider(names_from=sample, values_from=expression) %>%
    rename_with(function(x){if_else(x=="gene", true=x, false=paste0("NA", str_sub(x, 1, 5)))})
  expression.mat <- expression %>%
    column_to_rownames("gene") %>%
    as.matrix
  
  # expression normalization
  gene_names <- rownames(expression.mat)
  inds <- colnames(expression.mat)
  # normalize each gene to have zero mean, unit variance across individuals
  expression.mat <- apply(expression.mat, 1, invnorm_transform) %>% t
  # re-insert gene names and convert to tibble
  rownames(expression.mat) <- gene_keepers
  colnames(expression.mat) <- inds
  expression <- as_tibble(expression.mat, rownames="gene")
  
  # expression PCs - produces a scree plot and a PC biplot, and saves the covariates
  pcomp <- princomp(expression.mat)
  
  # viz PC biplot, scree, and num cells for the cell type
  png(file=paste0(plot_prefix, "/", d, "/qc.png"), width=480*3, height=480)
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
    write_tsv(paste0(table_prefix, "/", d, "/all_pcs.tsv"))
}