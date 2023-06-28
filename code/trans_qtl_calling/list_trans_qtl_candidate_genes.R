library(nullranges)
library(vroom)
library(tidyverse)
set.seed(1234)

# tests_list <- "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt"
# af_loc <- "data/genotypes/af_all.frq"
# eb_sighits_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.tsv"
# eb_sighits_overlap_loc <- "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.full_gtex_overlap.bed"
# gtf_loc <- "data/gencode/gencode.hg38.filtered.gtf"
# gmt_file <- "data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
# 
# go_geneset <- "GOBP_TISSUE_DEVELOPMENT"
  
gtf_loc <- snakemake@input[['gtf']]
gmt_file <- snakemake@input[['gmt']]
go_geneset <- snakemake@wildcards[['gs']]
output_loc <- snakemake@output[['gene_set']]

## List the genes in the gene set of interest
gmt_lines <- readLines(gmt_file)
gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, "\t")))
gmt_list2 <- lapply(gmt_list, function(x) list(x[1], x[2], x[3:length(x)]))
gmt_df <- as_tibble(do.call("rbind", gmt_list2)) %>%
  mutate(across(V1:V2, unlist)) %>%
  dplyr::rename(geneset=V1, link=V2, genes=V3) %>%
  filter(geneset==go_geneset) %>%
  select(-c(link)) %>%
  unnest(genes)

## Get a mapping from HGNC to ENSG
gencode <- vroom(gtf_loc, col_select=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute", "hgnc", "ensg"))

gs_genes <- inner_join(select(gmt_df, genes),
                       select(gencode, c(hgnc, ensg)),
                       by=c("genes"="hgnc")) %>%
  select(ensg) %>%
  write_tsv(output_loc, col_names=F)
