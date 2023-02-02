library(tidyverse)
library(vroom)

pull_gene_type <- function(attr) {
  str_split(attr, "\"")[[1]][6]
}
pull_gene_name <- function(attr) {
  str_split(attr, "\"")[[1]][8]
}
pull_gene_id <- function(attr) {
  str_split(attr, "\"")[[1]][2]
}

protein_coding_genes <- vroom("/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf", 
                 col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=5) %>%
  filter(seqname %in% paste0("chr", seq(1, 22))) %>%
  filter(feature == "gene") %>%
  mutate(type=map_chr(attribute, pull_gene_type)) %>%
  mutate(ensg=map_chr(attribute, pull_gene_id)) %>%
  mutate(hgnc=map_chr(attribute, pull_gene_name)) %>%
  filter(type=="protein_coding") %>%
  select(hgnc) %>%
  rename(gene=hgnc) %>%
  write_tsv("/project2/gilad/jpopp/ebQTL/data/fca/protein_coding_genes.tsv")