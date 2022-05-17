### 
# Aggregate pseudobulk within each cell type
# 
# Annotation prefix allows for multiple approaches to cell type annotation to be efficiently used in the workflow (ie multiple resolution clusterings)
# Uses TMM-normalization to account for (individual, cell type) samples having different depths
# Removes any samples with less than 10 (or `min_cells_per_sample` value) cells
###

library(anndata)
library(tidyverse)
library(Matrix)
library(Matrix.utils)
library(edgeR)
set.seed(2021)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")
#use_condaenv("/project2/gilad/jpopp/ebQTL/.snakemake/conda/980130bc9d8d3f8069c208a84393fb39")

# Read in the input files from Snakefile
adata_loc <- snakemake@input[["anndata"]]
pseudobulk_loc <- snakemake@input[["pseudobulk"]]
counts_loc <- snakemake@input[["counts"]]
sample_summary_loc <- as.character(snakemake@output[["sample_summary"]])
celltype_summary_loc <- as.character(snakemake@output[["celltype_summary"]])

# Read in the cell types and annotation information from Snakefile
extract_celltype <- function(s) {
  str_split(s, "/")[[1]] %>% rev %>% magrittr::extract2(2)
} 
select_celltype_files <- function(l) {
  l[sapply(l, function(s){!grepl("sample", s, fixed=TRUE)})]
}
cell.types <- snakemake@output %>% 
  select_celltype_files %>% 
  map_chr(extract_celltype) %>%
  unique
annotation_prefix <- snakemake@output %>%
  tail(1) %>%
  str_match(".*/") %>%
  as.character
print(paste0("here is the snakemake output: \n", snakemake@output))
print(paste0("\n \n here are the cell types: \n", cell.types))

# Set cutoff
min_cells_per_sample <- 5

# HELPER FUNCTIONS
counts2cpm <- function(c, lognorm=F) {
  g <- c$gene
  c <- select(c, !gene)
  head(c, n=5)
  cpm <- DGEList(counts=c) %>% calcNormFactors(method="TMM") %>% cpm(log=lognorm) %>%
    as_tibble %>% mutate(gene=g, .before=1) %>% arrange(gene)
  cpm
}

# # PSEUDOBULK
# adata <- anndata::read_h5ad(adata_loc, backed=NULL)
# counts <- adata$X
# 
# # aggregate pseudobulk by cell type
# type.ind <- paste(adata$obs[['individual']], adata$obs[['type']], sep="_")
# counts.type <- counts %>% aggregate.Matrix(type.ind, fun="sum") %>% t %>% 
#   as.matrix %>% as_tibble(rownames="gene") %>% arrange(gene)
counts.type <- read_tsv(pseudobulk_loc) %>% 
  arrange(gene)
type.ind <- colnames(counts.type)[-c(1)]

# create a summary file with the number of cells, number of reads per (individual, cell type) sample
# cell_counts <- as_tibble(table(type.ind)) %>%
#   rename(ind_type=type.ind, n_cells=n) %>%
#   separate(ind_type, into=c("individual", "type"), remove=F)
cell_counts <- read_tsv(counts_loc) %>%
    separate(ind_type, into=c("individual", "type"), sep="_", remove=F)

summary <- counts.type %>%
  column_to_rownames("gene") %>% 
  colSums %>% as_tibble(rownames="ind_type") %>%
  dplyr::rename(n_umi=value) %>%
  right_join(cell_counts, by="ind_type") %>%
  arrange(n_cells) %>%
  rowwise() %>%
  mutate(dropped=if_else(n_cells < min_cells_per_sample, "T", "F")) %>%
  write_tsv(sample_summary_loc)

# create a second summary file with the number of samples present per cell type
n_samples <- summary %>%
  filter(dropped=="F") %>%
  dplyr::count(type) %>%
  write_tsv(celltype_summary_loc)

# subset the expression data to the samples with at least 10 [see threshold above] cells
samples.subset <- summary %>%
  filter(dropped=="F") %>%
  arrange(ind_type)
counts.type <- counts.type %>%
  select(gene, all_of(samples.subset$ind_type))

# normalize, compute expression pcs, and save important information
# NA prepended to sample names in accordance with genotype data
for (d in cell.types) {
  print(d)
  expression <- counts.type %>% 
    select(c(gene, ends_with(paste0("_", d)))) %>% 
    rename_with(function(x){if_else(x=="gene", true=x, false=paste0("NA", str_sub(x, 1, 5)))}) %>% 
    counts2cpm(lognorm=T) %>% 
    write_tsv(paste0(annotation_prefix, d, "/expression.tsv"))
  covariates <- regular.pca(expression, 5)$u %>% 
    column_to_rownames("sample") %>% t %>% as_tibble(rownames="covariate") %>%
    write_tsv(paste0(annotation_prefix, d, "/covariates.tsv"))
  individuals <- samples.subset %>%
    filter(type==!!d) %>%
    select(individual) %>%
    mutate(individual=paste0("NA", individual)) %>%
    write_tsv(paste0(annotation_prefix, d, "/individuals.tsv"), col_names=FALSE)
}