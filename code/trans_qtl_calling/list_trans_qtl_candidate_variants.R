library(nullranges)
library(vroom)
library(tidyverse)
set.seed(1234)

eb_sighits_loc <- snakemake@input[['eb_hits']]
eb_sighits_overlap_loc <- snakemake@input[['eb_gtex_overlap']]
gmt_file <- snakemake@input[['gmt']]
bim_file <- snakemake@input[['bim']]

go_geneset <- snakemake@wildcards[['gs']]

variant_candidates_info_loc <- snakemake@output[['candidate_info']]
variant_candidates_loc <- snakemake@output[['candidates']]
matchers_loc <- snakemake@output[['match_details']]

## List the eQTLs for genes with no GTEx overlap that belong to this gene set
gtex_overlap <- vroom(eb_sighits_overlap_loc, col_names=c("CHR", "START", "STOP", "EB_ENSG", "EB_HGNC", 
                                                          "RSID", "CELLTYPE", "CHR1", "START1", 
                                                          "STOP1", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT"))
eb_all_hits <- vroom(eb_sighits_loc)

overlap_egenes <- unique(gtex_overlap$EB_HGNC)
all_egenes <- unique(eb_all_hits$phenotype_id)
novel_egenes <- setdiff(all_egenes, overlap_egenes)

## List the genes in the gene set of interest with and without GTEx overlap
gmt_lines <- readLines(gmt_file)
gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, "\t")))
gmt_list2 <- lapply(gmt_list, function(x) list(x[1], x[2], x[3:length(x)]))
gmt_df <- as_tibble(do.call("rbind", gmt_list2)) %>%
  mutate(across(V1:V2, unlist)) %>%
  dplyr::rename(geneset=V1, link=V2, genes=V3)
gmt_mat <- gmt_df %>%
  select(-c(link)) %>%
  unnest(genes) %>%
  filter(genes %in% all_egenes) %>%
  pivot_wider(names_from = genes, values_from = genes, values_fn = length, values_fill = 0) %>%
  column_to_rownames("geneset") %>% 
  as.matrix %>% t

gs_genes <- rownames(gmt_mat)[which(gmt_mat[,go_geneset]==1)]
novel_egenes_gs <- intersect(gs_genes, novel_egenes)

## List the eQTLs associated with the 'novel' eGenes 
pathway_novel <- filter(eb_all_hits, phenotype_id %in% novel_egenes_gs) 

# Load the bim file to make the variant format compatible w GTEx
bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS_BP", "ALLELE_1", "ALLELE_2", "variant_id"))
pathway_novel_qtls <- left_join(pathway_novel, bim, by=c("variant_id")) %>%
  mutate(variant_id_gtex = paste(`#CHR`, POS_BP, ALLELE_2, ALLELE_1, "b38", sep="_"))

# Save just this set of variants 
pathway_novel_qtls %>%
  dplyr::select(variant_id_gtex, phenotype_id) %>%
  write_tsv(variant_candidates_info_loc)

pathway_novel_qtls %>%
  dplyr::select(variant_id) %>%
  write_tsv(variant_candidates_loc, col_names=F)
