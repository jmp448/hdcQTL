library(tidyverse)
library(mashr)
library(vroom)
library(TreeBH)
set.seed(1234)

# Get input & output file names from snakemake
# cell.types <- read_tsv("data/benchmark_specificity_methods/eb_cellid/eb_cellid.pseudobulk_tmm.tsv", n_max=0, col_select=-c(gene)) %>%
#   colnames %>%
#   sapply(function(s){str_sub(s, start=7)}) %>%
#   unique()
# qtl_files <- paste0("results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/basic/", cell.types, "/8pcs/matrixeqtl.cis_qtl_pairs.all.tsv")

qtl_files <- snakemake@input
hierarchy <- as.character(snakemake@wildcards[['hierarchy']])
qtl_master_loc <- as.character(snakemake@output[['all_qtls']])
qtl_significance_loc <- as.character(snakemake@output[['qtl_significance']])

# helper functions
pull_type <- function(s) {
  str_split(s, "/")[[1]][6]
}

# compile a list of QTL test p-values across all tests
qtls <- vroom(qtl_files, id="path", col_select=c(SNP, gene, `p-value`, path)) %>%
  mutate(context=sapply(path, pull_type), .keep="unused") %>%
  unite(context_gene, c(context, gene), sep="--", remove=FALSE) %>%
  unite(test, c(context, gene, SNP), sep="--", remove=FALSE) %>%
  write_tsv(qtl_master_loc)

# create dictionaries mapping each context, gene, and test to an integer 
if (hierarchy == "context-gene-test") {
  context_dict <- tibble(context=unique(qtls$context)) %>%
    rowid_to_column(var="context_id")
  context_gene_dict <- tibble(context_gene=unique(qtls$context_gene)) %>%
    rowid_to_column(var="context_gene_id")
  test_dict <- tibble(test=unique(qtls$test)) %>%
    rowid_to_column(var="test_id")
  stopifnot(nrow(test_dict)==nrow(qtls)) # last column must be 1:N
} else if (hierarchy == "gene-context-test") {
  gene_dict <- tibble(gene=unique(qtls$gene)) %>%
    rowid_to_column(var="gene_id")
  context_gene_dict <- tibble(context_gene=unique(qtls$context_gene)) %>%
    rowid_to_column(var="context_gene_id")
  test_dict <- tibble(test=unique(qtls$test)) %>%
    rowid_to_column(var="test_id")
}

if (hierarchy == "context-gene-test") {
  qtls_indexed <- qtls %>%
    left_join(context_dict, by="context") %>%
    left_join(context_gene_dict, by="context_gene") %>%
    left_join(test_dict, by="test") 
  
  treebh_all <- qtls_indexed %>%
    select(c(context_id, context_gene_id, test_id, `p-value`))%>%
    as.matrix
} else if (hierarchy == "gene-context-test") {
  qtls_indexed <- qtls %>%
    left_join(gene_dict, by="gene") %>%
    left_join(context_gene_dict, by="context_gene") %>%
    left_join(test_dict, by="test") 
  
  treebh_all <- qtls_indexed %>%
    select(c(gene_id, context_gene_id, test_id, `p-value`)) %>%
    as.matrix
}

# Run TreeBH
treebh_selections <- get_TreeBH_selections(pvals=treebh_all[,4], groups=treebh_all[,1:3], q=c(0.05, 0.05, 0.05))

# Add significance testing results back to where we're keeping track of all tests
if (hierarchy == "context-gene-test") {
  qtls_indexed <- qtls_indexed %>%
    left_join(unique(tibble(context_id=treebh_all[,1], context_significance=as.logical(treebh_selections[,1]))), by="context_id") %>%
    left_join(unique(tibble(context_gene_id=treebh_all[,2], context_gene_significance=as.logical(treebh_selections[,2]))), by="context_gene_id") %>%
    left_join(unique(tibble(test_id=treebh_all[,3], test_significance=as.logical(treebh_selections[,3]))), by="test_id")
} else if (hierarchy == "gene-context-test") {
  qtls_indexed <- qtls_indexed %>%
    left_join(unique(tibble(gene_id=treebh_all[,1], gene_significance=as.logical(treebh_selections[,1]))), by="gene_id") %>%
    left_join(unique(tibble(context_gene_id=treebh_all[,2], context_gene_significance=as.logical(treebh_selections[,2]))), by="context_gene_id") %>%
    left_join(unique(tibble(test_id=treebh_all[,3], test_significance=as.logical(treebh_selections[,3]))), by="test_id")
}

qtl_significance <- qtls_indexed %>%
  select(-ends_with("_id")) %>%
  write_tsv(qtl_significance_loc)
