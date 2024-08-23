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

d_stats_loc <- as.character(snakemake@output[["dstats"]])
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
  dplyr::select(-c(gene)) %>%
  colnames %>%
  sort
pseudobulk <- relocate(pseudobulk, all_of(c("gene", sample_order)))

# Identify low-quality samples & sort sample summary
samples_subset <- read_tsv(sample_summary_loc) %>%
  arrange(ind_type) %>%
  filter(!dropped)

pseudobulk <- pseudobulk %>%
  dplyr::select(c(gene, intersect(colnames(.), samples_subset$ind_type)))

d_statistics <- tibble(sample=character(), D=numeric())
# Prep for QTL calling
for (d in c(cell.types)) {
  # get expression
  expression <- pseudobulk %>% 
    dplyr::select(c(gene, ends_with(paste0("_", d)))) %>% 
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
  
  # expression normalization
  # 1. TMM normalization to get log TMM-normalized counts
  expression.mat <- DGEList(counts=expression.mat) %>% calcNormFactors(method="TMM") %>% cpm(log=T)
  # 2. re-insert gene names and convert to tibble
  rownames(expression.mat) <- gene_keepers
  colnames(expression.mat) <- inds
  
  sample_corrs <- cor(expression.mat)
  diag(sample_corrs) <- NA # don't include self-correlation [always 1] for d-statistics
  d_stats_ct <- rowMeans(sample_corrs, na.rm=T)
  d_statistics <- bind_rows(d_statistics, tibble(sample=paste0(names(d_stats_ct), "_", d), D=d_stats_ct))
}

write_tsv(d_statistics, d_stats_loc)

# 
d_statistics <- read_tsv("/project2/gilad/jpopp/ebQTL/data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/d_statistics.tsv")
sample_summary <- read_tsv("data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/sample_summary_manual.tsv") %>% 
  mutate(ind_type=paste0("NA", ind_type)) %>%
  left_join(d_statistics, by=c("ind_type"="sample")) %>%
  arrange(type, desc(n_cells_unfiltered)) 

