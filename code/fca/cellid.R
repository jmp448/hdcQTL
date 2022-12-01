library(SingleCellExperiment)
library(CelliD)
library(tidyverse)
library(vroom)

# Load subsampled fetal cell atlas data
fca_counts <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/counts.sampled.rds")
fca_cells <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/cell_metadata.rds")
fca_genes <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/gene_metadata.rds") %>%
  mutate_all(as.character)

# Subset cells to those in the subsampled data
fca_cells <- fca_cells %>%
  filter(sample %in% colnames(fca_counts))

# Subset genes to protein-coding genes without duplicated HGNC entries
pull_gene_type <- function(attr) {
  str_split(attr, "\"")[[1]][6]
}
pull_gene_name <- function(attr) {
  str_split(attr, "\"")[[1]][8]
}
pull_gene_id <- function(attr) {
  str_split(attr, "\"")[[1]][2]
}

gencode <- vroom("/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf", 
                 col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=5) %>%
  filter(seqname %in% paste0("chr", seq(1, 22))) %>%
  filter(feature == "gene") %>%
  mutate(type=map_chr(attribute, pull_gene_type)) %>%
  mutate(ensg=map_chr(attribute, pull_gene_id)) %>%
  mutate(hgnc=map_chr(attribute, pull_gene_name)) %>%
  filter(type=="protein_coding")

fca_genes <- mutate(fca_genes, gene_id_short=str_extract(gene_id, "[^.]+"))
fca_pc_genes <- filter(fca_genes, 
                       (gene_short_name %in% gencode$hgnc) & (gene_id_short %in% gencode$ensg))
hgnc_duplicates <- fca_pc_genes %>%
  group_by(gene_short_name) %>%
  count() %>%
  filter(n > 1) %>%
  pull(gene_short_name)
fca_pc_genes <- fca_pc_genes %>%
  filter(!gene_short_name %in% hgnc_duplicates)

# Make the counts matrix match the metadata
fca_counts <- fca_counts[fca_pc_genes$gene_id,fca_cells$sample]
rownames(fca_counts) <- fca_pc_genes$gene_short_name

# Create SingleCellExperiment object
fca <- SingleCellExperiment(list(counts=fca_counts),
                            colData=fca_cells)
saveRDS(fca, "/project2/gilad/jpopp/ebQTL/data/fca/counts.sampled.sce")
