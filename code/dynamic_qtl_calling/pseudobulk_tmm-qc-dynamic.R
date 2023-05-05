### 
# Prep for eQTL calling by measuring expression covariates and cell line PCs
# 
# Annotation prefix allows for multiple approaches to cell type annotation to be efficiently used in the workflow (ie multiple resolution clusterings)
# Relies on having previously done aggregation in a notebook like `analysis/highpass/aggregation_*`
#
###

library(tidyverse)
library(Matrix)
library(matrixStats)
library(edgeR)
set.seed(2021)
source("/project2/gilad/jpopp/sc-dynamic-eqtl/code/cell_line_pca.R")

pseudobulk_loc <- snakemake@input[["pseudobulk"]]
sample_summary_loc <- as.character(snakemake@input[["sample_summary"]])
sample_summary_manual_loc <- as.character(snakemake@output[["sample_summary_manual"]])
table_prefix <- snakemake@params[['table_prefix']]
plot_prefix <- snakemake@params[['fig_prefix']]

# Inverse normal function
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

pseudobulk <- pseudobulk %>%
  select(c(gene, intersect(colnames(.), samples_subset$ind_type)))

# Set up sample summary for manual check
if (!file.exists(sample_summary_manual_loc)) {
  file.copy(sample_summary_loc, sample_summary_manual_loc)
}

# Set up directories if they don't exist
if(!dir.exists(plot_prefix)) {
  dir.create(plot_prefix, recursive=T)
  print("Plot directory created.")
} else {
  print("Plot directory already existed.")
}

if(!dir.exists(table_prefix)) {
  dir.create(table_prefix, recursive=T)
  print("Table directory created.")
} else {
  print("Table directory already existed.")
}

expression.mat <- pseudobulk %>%
  column_to_rownames("gene") %>%
  as.matrix

# gene filter - nonzero variance, nonzero median expression
gene_vars <- rowVars(expression.mat) > 0
gene_nonzeros <- rowMedians(expression.mat) > 10 
gene_keepers <- rownames(expression.mat)[gene_vars & gene_nonzeros]
samples <- colnames(expression.mat)
expression.mat <- expression.mat[gene_keepers,]

# expression normalization
# 1. TMM normalization to get log TMM-normalized counts
expression.mat <- DGEList(counts=expression.mat) %>% calcNormFactors(method="TMM") %>% cpm(log=T)
# 2. normalize each gene to have zero mean, unit variance across individuals
expression.mat <- apply(expression.mat, 1, invnorm_transform) %>% t
# 3. re-insert gene names and convert to tibble
rownames(expression.mat) <- gene_keepers
colnames(expression.mat) <- samples
expression <- as_tibble(expression.mat, rownames="gene") 

# Expression PCs - sanity check / sample QC
pcomp <- regular.pca(expression)

expression_pcs <- pcomp$u %>% 
  column_to_rownames("sample") %>% t %>% as_tibble(rownames="covariate") %>%
  write_tsv(paste0(table_prefix,"all_expression_pcs.tsv"))

# visualize PC biplot and scree
png(file=paste0(plot_prefix, "expression_pc_qc.png"), width=480*3, height=480)
par(mfrow=c(1,2))
plot(pcomp$d ^ 2 / sum(pcomp$d ^ 2), xlab="PC", ylab="Variance Explained")
plot(pcomp$u$PC1, pcomp$u$PC2, 
     col=colorRampPalette(c("gray", "red"))(nrow(samples_subset))[rank(samples_subset$pseudotime)], 
     xlab="PC1", ylab="PC2")
dev.off()

# Cell line PCs - covariates for QTL calling
cell_line_pcs <- cell.line.pca(expression,npc=15)
cell_line_pcs$cell.line.pcs %>% 
  write_tsv(paste0(table_prefix,'all_cell_line_pcs.tsv'))

# visualize cell line PC biplot and scree
png(file=paste0(plot_prefix, "cell_line_pc_qc.png"), width=480*3, height=480)
par(mfrow=c(1,2))
plot(cell_line_pcs$pve ^ 2 / sum(cell_line_pcs$pve ^ 2), xlab="PC", ylab="Variance Explained")
plot(cell_line_pcs$cell.line.pcs$PC_1,cell_line_pcs$cell.line.pcs$PC_2)
dev.off()
