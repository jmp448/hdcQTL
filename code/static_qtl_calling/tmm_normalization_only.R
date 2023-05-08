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

# Read in the input files from Snakefile
pseudobulk_loc <- "data/static_qtl_calling/elorbany_cmstages/pseudobulk_tmm/elorbany_cmstages.pseudobulk_tmm.tsv"
sample_summary_loc="data/static_qtl_calling/elorbany_cmstages/pseudobulk_tmm/sample_summary.tsv"
celltypes_loc="data/static_qtl_calling/elorbany_cmstages/pseudobulk_tmm/elorbany_cmstages.pseudobulk_tmm.tsv"

tmm_expression_loc <- "temp/elorbany_cmstages.pseudobulk_tmm.tsv"

pseudobulk_loc <- snakemake@input[['pseudobulk']]
sample_summary_loc <- snakemake@input[['sample_summary_manual']]
celltypes_loc <- snakemake@input[['celltypes']]

full_expression_loc <- snakemake@output[['all']]
median_expression_loc <- snakemake@output[['median']]

# Get a list of cell types
cell.types <- read_tsv(celltypes_loc, n_max=0, col_select=-c(1)) %>%
  colnames() %>%
  sapply(function(s){str_split(s, "_")[[1]][2]}) %>% 
  unique()

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

pseudobulk <- pseudobulk %>%
  select(c(gene, intersect(colnames(.), samples_subset$ind_type)))

for (d in c(cell.types)) {
  expression <- pseudobulk %>% 
    select(c(gene, ends_with(paste0("_", d))))
  expression.mat <- expression %>%
    column_to_rownames("gene") %>%
    as.matrix
  
  # TMM normalization to get log TMM-normalized counts
  expression.mat <- DGEList(counts=expression.mat) %>% calcNormFactors(method="TMM") %>% cpm(log=F)
  expression <- as_tibble(expression.mat, rownames="gene")
  
  if (d == cell.types[[1]]) {
    full_expression <- expression
  } else {
    full_expression <- full_join(full_expression, expression, by="gene")
  }
}

full_expression %>%
  write_tsv(full_expression_loc)

median_expression <- full_expression %>%
  pivot_longer(!gene, names_to="sample", values_to="expression") %>%
  separate(sample, into=c("donor", "stage"), sep="_") %>%
  group_by(stage, gene) %>%
  summarize(median_expression=median(expression)) %>%
  pivot_wider(id_cols=gene, names_from=stage, values_from=median_expression) %>%
  write_tsv(median_expression_loc)