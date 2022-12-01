library(zellkonverter)
library(scater)
library(scran)
library(Matrix.utils)
library(edgeR)
library(tidyverse)

adata <- readH5AD("/project2/gilad/jpopp/ebQTL/data/single_cell_objects/highpass/cellid_annotated.ipscfilt_only.h5ad")

sce <- adata

# Check whether gene filtering is necessary
counts_cpm <- cpm(counts(adata), log = TRUE)
mean_logcpm <- rowMeans(counts_cpm)
summary(mean_logcpm)

# Perform scran normalization
clusts <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusts)

# Get normalized expression
sce <- logNormCounts(sce)

# Get average expression of each gene in each individual
exp <- sce@assays@data@listData$logcounts
inds <- sce@colData$donor_id_x
pseudobulk <- as_tibble(as.matrix(t(aggregate.Matrix(t(exp), groupings=inds, fun="mean"))), rownames="gene")

# To mirror our other pseudobulk objects, I want the rework the column (individual) names
rename_inds <- function(i) {
  paste0(str_sub(i, 3), "_IPSC")
}
pseudobulk <- pseudobulk %>% rename_with(rename_inds, .cols=!gene)

# save to a TSV
write_tsv(pseudobulk, "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/ebqtl_ipscfilt.pseudobulk-scran.tsv")

# write a sample summary - ind_type, n_umi, individual, type, n_cells, dropped
sample_summary <- as.data.frame(sce@colData@listData) %>%
  select(donor_id_x) %>%
  mutate(ind_type=rename_inds(donor_id_x), .keep="unused") %>%
  count(ind_type, name="n_cells") %>%
  mutate(dropped=n_cells < 5) %>%
  separate(ind_type, into=c("individual", "type"), sep="_", remove=F) %>%
  mutate(n_umi=NA) %>%
  relocate(ind_type, n_umi, individual, type, n_cells, dropped) %>%
  write_tsv("/project2/gilad/jpopp/ebQTL/data/static/ebqtl_ipscfilt/pseudobulk-scran/sample_summary.tsv")
