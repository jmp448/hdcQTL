library(tidyverse)
library(vroom)

# gmt_file <- "/project2/gilad/jpopp/ebQTL/data/gene_sets/c5.go.bp.v2023.1.Hs.symbols.gmt"
gmt_file <- snakemake@input[['gmt']]
gtf_file <- snakemake@input[['gtf']]
gs <- snakemake@wildcards[['geneset']]
gs_loc <- snakemake@output[['gs']]

gene_dict <- vroom(gtf_file)
gmt_lines <- readLines(gmt_file)

gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, "\t")))
gmt_list2 <- lapply(gmt_list, function(x) list(x[1], x[2], x[3:length(x)]))

gmt_df <- as_tibble(do.call("rbind", gmt_list2)) %>%
  mutate(across(V1:V2, unlist)) %>%
  dplyr::rename(geneset=V1, link=V2, genes=V3)

gs_list <- tibble("hgnc"=unlist(filter(gmt_df, geneset == str_replace_all(gs, "-", "_"))$genes)) %>%
  left_join(select(gene_dict, c(hgnc, ensg)), by="hgnc") %>%
  write_tsv(gs_loc)

# Manually combine hepaticobiliary system and kidney
# hep <- tibble("hgnc"=unlist(filter(gmt_df, geneset == "GOBP_HEPATICOBILIARY_SYSTEM_DEVELOPMENT")$genes)) %>%
#   left_join(select(gene_dict, c(hgnc, ensg)), by="hgnc")
# renal <- tibble("hgnc"=unlist(filter(gmt_df, geneset == "GOBP_RENAL_SYSTEM_DEVELOPMENT")$genes)) %>%
#   left_join(select(gene_dict, c(hgnc, ensg)), by="hgnc")
# hep_and_renal <- bind_rows(hep, renal) %>%
#   distinct() %>%
#   write_tsv("data/ldsc/genesets/GOBP-HEPATICOBILIARY-AND-RENAL-SYSTEM-DEVELOPMENT.tsv")
