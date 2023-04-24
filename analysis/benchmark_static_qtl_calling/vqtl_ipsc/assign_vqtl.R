library(scater)
library(scran)
library(Matrix.utils)

eset <- readRDS("/project2/gilad/jpopp/ebQTL/data/vqtl_paper/GSE118723_eset.rds")
dim(eset)

eset <- eset[fData(eset)$source %in% c("H. sapiens") , ]
dim(eset) 

# load quality cells, quality genes
quality_cells <- read_tsv("/project2/gilad/jpopp/ebQTL/data/vqtl_paper/GSE118723_quality-single-cells.txt",
                          col_names=c("cell", "qual"))
cell_keepers <- quality_cells %>%
  filter(qual==T) %>%
  pull(cell)

quality_genes <- read_tsv("/project2/gilad/jpopp/ebQTL/data/vqtl_paper/GSE118723_genes-pass-filter.txt",
                          col_names=c("gene", "qual"))
gene_keepers <- quality_genes %>%
  filter(qual==T) %>%
  pull(gene)

# Get the raw counts
counts <- eset@assayData$exprs

# Get the metadata
metadata <- eset@phenoData@data

# Create a single cell experiment object
sce <- SingleCellExperiment(list(counts=counts), colData=metadata)

# Compute library size factors
lib.sf.sce <- librarySizeFactors(sce)

# Get normalized expression
sce <- logNormCounts(sce, size_factors = lib.sf.sce)

# Filter to good genes & cells
sce <- sce[gene_keepers, cell_keepers]

saveRDS(sce, "/project2/gilad/jpopp/ebQTL/data/benchmark_static_qtl_calling/vqtl_ipsc/vqtl_ipsc.sce")

# Get average & summed expression of each gene in each individual
exp_norm <- sce@assays@data@listData$logcounts
exp_raw <- sce@assays@data@listData$counts
inds <- sce@colData$chip_id
pseudobulk_scran <- as_tibble(as.matrix(t(aggregate.Matrix(t(exp_norm), groupings=inds, fun="mean"))), rownames="gene")
pseudobulk_tmm <- as_tibble(as.matrix(t(aggregate.Matrix(t(exp_raw), groupings=inds, fun="sum"))), rownames="gene")

# To mirror our other pseudobulk objects, I want the rework the column (individual) names
rename_inds <- function(i) {
  paste0(str_sub(i, 3), "_IPSC")
}
pseudobulk_scran <- pseudobulk_scran %>% rename_with(rename_inds, .cols=!gene)
pseudobulk_tmm <- pseudobulk_tmm %>% rename_with(rename_inds, .cols=!gene)

# We also want to replace the ENSG gene names with HGNC gene names
pull_ensg <- function(attr) {
  str_split(attr, "\"")[[1]][3]
}

pull_hgnc <- function(attr) {
  str_split(attr, "\"")[[1]][15]
}

gencode <- read_tsv("/project2/gilad/jpopp/ebQTL/data/gencode/gencode.hg38.filtered.gtf") %>%
  mutate(ensg=map_chr(attribute, pull_ensg)) %>%
  mutate(hgnc=map_chr(attribute, pull_hgnc))

gene_dict <- dplyr::select(gencode, c(ensg, hgnc)) %>%
  filter(ensg %in% pseudobulk_scran$gene) 

pseudobulk_scran <- pseudobulk_scran %>%
  left_join(gene_dict, by=c("gene"="ensg")) %>%
  select(!gene) %>%
  dplyr::rename(gene=hgnc) %>%
  relocate(gene, .before=1) %>%
  drop_na

pseudobulk_tmm <- pseudobulk_tmm %>%
  left_join(gene_dict, by=c("gene"="ensg")) %>%
  select(!gene) %>%
  dplyr::rename(gene=hgnc) %>%
  relocate(gene, .before=1) %>%
  drop_na

# save to a TSV
write_tsv(pseudobulk_scran, "/project2/gilad/jpopp/ebQTL/data/benchmark_static_qtl_calling/vqtl_ipsc/pseudobulk_scran/vqtl_ipsc.pseudobulk_scran.tsv")
write_tsv(pseudobulk_tmm, "/project2/gilad/jpopp/ebQTL/data/benchmark_static_qtl_calling/vqtl_ipsc/pseudobulk_tmm/vqtl_ipsc.pseudobulk_tmm.tsv")

# write a sample summary - ind_type, n_umi, individual, type, n_cells, dropped
sample_summary <- as.data.frame(sce@colData@listData) %>%
  select(chip_id) %>%
  mutate(ind_type=rename_inds(chip_id), .keep="unused") %>%
  dplyr::count(ind_type, name="n_cells_filtered") %>%
  mutate(dropped=n_cells_filtered < 5) %>%
  separate(ind_type, into=c("individual", "type"), sep="_", remove=F) %>%
  mutate(n_umi=NA) %>%
  relocate(ind_type, n_umi, individual, type, n_cells_filtered, dropped) %>%
  write_tsv("data/benchmark_static_qtl_calling/vqtl_ipsc/pseudobulk_scran/sample_summary.tsv")
  
sample_summary <- as.data.frame(sce@colData@listData) %>%
  select(chip_id) %>%
  mutate(ind_type=rename_inds(chip_id), .keep="unused") %>%
  dplyr::count(ind_type, name="n_cells_filtered") %>%
  mutate(dropped=n_cells_filtered < 5) %>%
  separate(ind_type, into=c("individual", "type"), sep="_", remove=F) %>%
  mutate(n_umi=NA) %>%
  relocate(ind_type, n_umi, individual, type, n_cells_filtered, dropped) %>%
  write_tsv("data/benchmark_static_qtl_calling/vqtl_ipsc/pseudobulk_tmm/sample_summary.tsv")