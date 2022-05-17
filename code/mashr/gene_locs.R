library(vroom)
library(tidyverse)

gtf_loc <- snakemake@input[['gtf_loc']]
filtered_gtf_loc <- snakemake@output[['gtf_loc']]
tss_loc <- snakemake@output[['tss_loc']]
bed_loc <- snakemake@output[['bed_loc']]

pull_gene_type <- function(attr) {
  str_split(attr, "\"")[[1]][6]
}

pull_gene_name <- function(attr) {
  str_split(attr, "\"")[[1]][8]
}

# load gencode data
gencode <- vroom(gtf_loc, col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=5) %>%
  filter(seqname %in% paste0("chr", seq(1, 22))) %>%
  filter(feature == "gene") %>%
  mutate(type=map_chr(attribute, pull_gene_type)) %>%
  mutate(hgnc=map_chr(attribute, pull_gene_name)) %>%
  filter(type=="protein_coding")

# 4 HGNC symbols are duplicated here - these have overlapping loci, and will be removed from our analysis
# TBCE, MATR3, HSPA14, GGT1
hgnc.dup <- filter(gencode, hgnc %in% names(table(gencode$hgnc)[table(gencode$hgnc)>1]))
gencode <- gencode %>% filter(!hgnc %in% hgnc.dup$hgnc)

# save a filtered GFF file
gtf.filt <- gencode %>%
  select(!c(type, hgnc)) %>%
  write_tsv(filtered_gtf_loc)

# save a tsv file of TSS for each gene with 1-indexing
tss.filt <- gencode %>%
  mutate(tss_start=if_else(strand=="+", start, end)) %>%
  select(hgnc, seqname, tss_start) %>%
  mutate(tss_end=tss_start)  %>%
  rename(chr=seqname) %>%
  write_tsv(tss_loc)

# save a bed file of TSS for each gene with [0,1) indexing
bed.filt <- gencode %>%
  mutate(tss=if_else(strand=="+", start, end)) %>%
  select(seqname, tss, hgnc) %>%
  mutate(bed.start=tss-1, .after=1) %>%
  write_tsv(bed_loc, col_names=F)
