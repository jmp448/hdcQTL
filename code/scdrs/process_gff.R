library(vroom)
library(tidyverse)

# filter the annotated genome to genes we found
gencode <- vroom(snakemake@input[[1]], col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=7) %>%
  filter(feature == "gene") %>%
  filter(seqname != "chrM") %>%
  dplyr::mutate(type=gsub(".*gene_type=([^;]+)[;].*", "\\1", attribute)) %>%
  dplyr::mutate(hgnc=gsub(".*gene_name=([^;]+)[;].*", "\\1", attribute)) %>%
  filter(type=="protein_coding")

# a few HGNC symbols are duplicated here - these have overlapping loci, and will be removed from our analysis
hgnc.dup <- filter(gencode, hgnc %in% names(table(gencode$hgnc)[table(gencode$hgnc)>1]))
gencode <- gencode %>% filter(!hgnc %in% hgnc.dup$hgnc)

# save a filtered GFF file
gff.filt <- gencode %>%
  select(!c(type, hgnc)) %>%
  write_tsv(snakemake@output[[1]])

# save a bed file of the gene locations with [0,1) indexing
bed.filt <- gencode %>%
  mutate(tss=if_else(strand=="+", start, end)) %>%
  mutate(tes=if_else(strand=="+", end, start)) %>%
  dplyr::select(seqname, tss, tes, hgnc) %>%
  mutate(tss=tss-1) %>%
  write_tsv(snakemake@output[[2]], col_names=F)

# save a bed file of TSS for each gene with [0,1) indexing
bed.tss.filt <- gencode %>%
  mutate(tss=if_else(strand=="+", start, end)) %>%
  dplyr::select(seqname, tss, hgnc) %>%
  mutate(bed.start=tss-1, .after=1) %>%
  write_tsv(snakemake@output[[3]], col_names=F)