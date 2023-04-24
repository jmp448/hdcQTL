library(tidyverse)
library(vroom)

gtf_loc <- snakemake@input[[1]]
pcg_list_loc <- snakemake@output[[1]]

# Some helper functions for reading the GTF file
pull_gene_type <- function(attr) {
  str_split(attr, "\"")[[1]][6]
}
pull_gene_name <- function(attr) {
  str_split(attr, "\"")[[1]][8]
}
pull_gene_id <- function(attr) {
  str_split(attr, "\"")[[1]][2]
}

# Read GTF file
protein_coding_genes <- vroom(gtf_loc, 
                 col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=5) %>%
  filter(seqname %in% paste0("chr", seq(1, 22))) %>%
  filter(feature == "gene") %>%
  mutate(type=map_chr(attribute, pull_gene_type)) %>%
  mutate(ensg=map_chr(attribute, pull_gene_id)) %>%
  mutate(hgnc=map_chr(attribute, pull_gene_name)) %>%
  filter(type=="protein_coding") %>%
  select(ensg, hgnc)

# Remove genes with duplicated HGNC
hgnc_duplicates <- protein_coding_genes %>%
  group_by(hgnc) %>%
  dplyr::count() %>%
  filter(n > 1) %>%
  pull(hgnc)
protein_coding_genes <- protein_coding_genes %>%
  filter(!hgnc %in% hgnc_duplicates) %>%
  write_tsv(pcg_list_loc)