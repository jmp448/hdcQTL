library(Matrix)
library(Matrix.utils)
library(scran)
library(SingleCellExperiment)
library(tidyverse)
library(CelliD)

set.seed(42)

fca_counts_loc="data/fca/counts.sampled.rds"
fca_cell_metadata_loc="data/fca/cell_metadata.rds"
fca_gene_metadata_loc="data/fca/gene_metadata.rds"
pcg_loc="data/gencode/gencode.hg38.filtered.gtf"

# Read snakemake inputs
fca_counts_loc <- snakemake@input[['counts']]
fca_cell_metadata_loc <- snakemake@input[['cell_metadata']]
fca_gene_metadata_loc <- snakemake@input[['gene_metadata']]
pcg_loc <- snakemake@input[['pc_genes']]

train_77_loc <- snakemake@output[['train_77']]
test_77_loc <- snakemake@output[['test_77']]
signatures_77_loc <- snakemake@output[['signatures_77']]

# Load subsampled fetal cell atlas data
fca_counts <- readRDS(fca_counts_loc)
fca_cells <- readRDS(fca_cell_metadata_loc)
fca_genes <- readRDS(fca_gene_metadata_loc) %>%
  mutate_all(as.character)

pc_genes <- read_tsv(pcg_loc)

# Subset cells to those in the subsampled data
fca_cells <- fca_cells %>%
  filter(sample %in% colnames(fca_counts))

# Subset genes to protein coding genes for normalization
rownames(fca_counts) <- str_extract(rownames(fca_counts), "[^.]+")
fca_pc_genes <- pc_genes %>%
  filter(ensg %in% rownames(fca_counts)) %>%
  arrange(ensg)

# Make the counts matrix match the metadata
fca_counts <- fca_counts[fca_pc_genes$ensg, fca_cells$sample] # subset genes & cells
rownames(fca_counts) <- fca_pc_genes$hgnc # switch genes to hgnc

# Create SingleCellExperiment object
fca <- SingleCellExperiment(list(counts=fca_counts),
                            colData=fca_cells)

# Revise cell type labels
fca$celltype <- str_extract(fca$Organ_cell_lineage, "[^-]+$")

# Update the dataset to contain even samples from these cell types, and a max of 5k cells per major cell type (regardless of organ)
celltype_subsetter <- as_tibble(colData(fca)) %>%
  dplyr::select(sample, celltype) %>%
  drop_na() %>%
  group_by(celltype) %>%
  add_count(celltype)
keepers_rare <- celltype_subsetter %>%
  filter(n <= 5000) %>%
  pull(sample)
keepers_sampled <- celltype_subsetter %>%
  filter(n > 5000) %>%
  slice_sample(n=5000) %>%
  pull(sample)
celltype_subsetter <- celltype_subsetter %>%
  filter(sample %in% c(keepers_rare, keepers_sampled)) 

fca_subsampled <- fca[,celltype_subsetter$sample]

# Now subset half the data for training, leaving half for testing
train_cells <- sample(colnames(fca_subsampled), size=round(ncol(fca_subsampled)/2), replace=F)
test_cells <- setdiff(colnames(fca_subsampled), train_cells)

fca_train <- fca_subsampled[,train_cells]
fca_test <- fca_subsampled[,test_cells]

# Normalize
libs <- librarySizeFactors(fca_train)
sizeFactors(fca_train) <- libs
fca_train <- logNormCounts(fca_train)

libs_test <- librarySizeFactors(fca_test)
sizeFactors(fca_test) <- libs_test
fca_test <- logNormCounts(fca_test)

fca_train.dec <- modelGeneVar(fca_train)
hvg.fca_train <- getTopHVGs(fca_train.dec, fdr.threshold=0.1)

# Learn signatures
fca_train <- RunMCA(fca_train, slot="logcounts", nmcs=50, features=hvg.fca_train)
fca_test <- RunMCA(fca_test, slot="logcounts", nmcs=50, features=hvg.fca_train)
fca_77_signatures <- GetGroupGeneSet(fca_train, group.by="celltype", n.features=100)

saveRDS(fca_train, train_77_loc)
saveRDS(fca_test, test_77_loc)
saveRDS(fca_77_signatures, signatures_77_loc)

