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
hvg_loc="data/fca/fca_subsampled_hvg.tsv"
pcg_loc="data/gencode/gencode.hg38.filtered.gtf"

# Read snakemake inputs
fca_counts_loc <- snakemake@input[['counts']]
fca_cell_metadata_loc <- snakemake@input[['cell_metadata']]
fca_gene_metadata_loc <- snakemake@input[['gene_metadata']]
pcg_loc <- snakemake@input[['pc_genes']]
hvg_loc <- snakemake@input[['hv_genes']]
fca_77_loc <- snakemake@output[['fca_77']]
signatures_77_loc <- snakemake@output[['signatures_77']]

# Load subsampled fetal cell atlas data
fca_counts <- readRDS(fca_counts_loc)
fca_cells <- readRDS(fca_cell_metadata_loc)
fca_genes <- readRDS(fca_gene_metadata_loc) %>%
  mutate_all(as.character)

hv_genes <- read_tsv(hvg_loc)
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
  select(sample, celltype) %>%
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

# Normalize
libs <- librarySizeFactors(fca_subsampled)
sizeFactors(fca_subsampled) <- libs
fca_subsampled <- logNormCounts(fca_subsampled)

# Learn signatures
fca_subsampled <- RunMCA(fca_subsampled, slot="logcounts", nmcs=50, features=hv_genes$gene)
fca_signatures <- GetGroupGeneSet(fca_subsampled, group.by="celltype", n.features=100)

saveRDS(fca_subsampled, fca_77_loc)
saveRDS(fca_signatures, signatures_77_loc)

