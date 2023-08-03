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
source("/project2/gilad/jpopp/ebQTL/code/dynamic_qtl_calling/clpca.R")

# sample_metadata_loc <- "/project2/gilad/katie/ebQTL/CombinedFormationAndCollectionMetadata_102andPilot_SWAPSANDCONTAMINATIONADDED_012522.csv"
# pseudobulk_loc <- "data/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/eb-cm_15binstrimmed.pseudobulk_tmm.tsv"
# sample_summary_loc <- "data/dynamic_qtl_calling/eb-cm_15binstrimmed/pseudobulk_tmm/sample_summary.tsv"

pseudobulk_loc <- snakemake@input[["pseudobulk"]]
sample_summary_loc <- as.character(snakemake@input[["sample_summary_manual"]])
sample_metadata_loc <- as.character(snakemake@input[["metadata"]])
table_prefix <- snakemake@params[['table_prefix']]
plot_prefix <- snakemake@params[['fig_prefix']]
pca_method <- as.character(snakemake@wildcards[['pca']])
raw_expression_loc <- snakemake@output[["raw_expression"]]
norm_expression_loc <- snakemake@output[["norm_expression"]]

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

# Prep for QTL calling
expression.mat <- pseudobulk %>%
  column_to_rownames("gene") %>%
  as.matrix

# gene filter - nonzero variance, nonzero median expression
gene_vars <- rowVars(expression.mat) > 0
gene_nonzeros <- rowMedians(expression.mat) > 10 
gene_keepers <- rownames(expression.mat)[gene_vars & gene_nonzeros]
samples <- colnames(expression.mat)
expression.mat <- expression.mat[gene_keepers,]

expression_raw <- as_tibble(expression.mat, rownames="gene") %>%
  write_tsv(raw_expression_loc)

# expression normalization
# 1. TMM normalization to get log TMM-normalized counts
expression.mat <- DGEList(counts=expression.mat) %>% calcNormFactors(method="TMM") %>% cpm(log=T)
# 2. normalize each gene to have zero mean, unit variance across individuals
expression.mat <- apply(expression.mat, 1, invnorm_transform) %>% t
# 3. re-insert gene names and convert to tibble
rownames(expression.mat) <- gene_keepers
colnames(expression.mat) <- samples
expression <- as_tibble(expression.mat, rownames="gene") %>%
  write_tsv(norm_expression_loc)

# expression PCs - produces a scree plot and a PC biplot, and saves the covariates
pcomp <- regular.pca(expression)

expression_pcs <- pcomp$u %>% 
  column_to_rownames("sample") %>% t %>% as_tibble(rownames="covariate") %>%
  write_tsv(paste0(table_prefix,"expression_pcs.tsv"))

# viz PC biplot and scree for all cells
png(file=paste0(plot_prefix, "expression_pc_qc_manual.png"), width=480*3, height=480)
par(mfrow=c(1,2))
plot(pcomp$d ^ 2 / sum(pcomp$d ^ 2), xlab="PC", ylab="Variance Explained")
plot(pcomp$u$PC1, pcomp$u$PC2, 
     col=colorRampPalette(c("gray", "red"))(nrow(samples_subset))[rank(samples_subset$pseudotime)], 
     xlab="PC1", ylab="PC2")
dev.off()

# Cell line PCs & covariates for QTL calling
metadata <- read_csv(sample_metadata_loc) %>%
  select(Line.True, sex) %>%
  mutate(ind=as.character(Line.True), .keep="unused") %>%
  distinct()
if (pca_method == "svd") {
  cell_line_pcs <- cell.line.pca(expression,npc=15)
} else {
  cell_line_pcs <- prob.cell.line.pca(expression,npc=15,pca_method)
}
covariates <- cell_line_pcs$cell.line.pcs %>%
  rename_with(str_replace, pattern="PC", replacement="CLPC") %>%
  left_join(metadata, by="ind") %>%
  relocate(sex, .after=ind) %>% 
  mutate(sex = if_else(sex=="F", 0, 1)) %>%
  write_tsv(paste0(table_prefix,'covariates.tsv'))

# Plot cell line PC biplot and scree
png(file=paste0(plot_prefix, "cell_line_pc_qc_manual.png"), width=480*3, height=480)
par(mfrow=c(1,2))
plot(cell_line_pcs$pve ^ 2 / sum(cell_line_pcs$pve ^ 2), xlab="PC", ylab="Variance Explained")
plot(cell_line_pcs$cell.line.pcs$PC_1,cell_line_pcs$cell.line.pcs$PC_2)
dev.off()

# Prepend NAs to sample names in accordance with genotype data
individuals <- samples_subset %>%
  mutate(individual=paste0("NA", individual), .keep="used") %>%
  write_tsv(paste0(table_prefix, "individuals.tsv"), col_names=FALSE)

# Save median pseudotime per sample for tensorqtl
pseudotime <- samples_subset %>%
  select(ind_type, pseudotime) %>%
  write_tsv(paste0(table_prefix, "pseudotime.tsv"))
