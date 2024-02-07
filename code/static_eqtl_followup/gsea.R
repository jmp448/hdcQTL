library(tidyverse)
library(Matrix.utils)
library(fgsea)
library(susieR)
library(pathways)
library(vroom)
library(gseasusie)

eb_bed <- snakemake@input[['eb_bed']]
gtex_overlap_bed <- snakemake@input[['gtex_bed']]
gmt_loc <- snakemake@input[['gmt']]

gsea_table_loc <- snakemake@output[['gsea_results']]

# Load eQTL sets - both the full HDC eQTL list, and the subset overlapping with GTEx
eb_distinct <- vroom(eb_bed) %>% 
  dplyr::select(c(EB_HGNC, EB_VARIANT_ID)) %>%
  distinct()

gtex_distinct <- vroom(gtex_overlap_bed, col_names=c("CHR", "START", "STOP", "EB_ENSG", "EB_HGNC", 
                                                           "RSID", "CELLTYPE", "CHR1", "START1", 
                                                           "STOP1", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT")) %>%
  dplyr::select(EB_HGNC, RSID) %>%
  distinct()

### Wrangle GO BP gene sets
gmt_lines <- readLines(gmt_loc)

gmt_list <- lapply(gmt_lines, function(x) unlist(strsplit(x, "\t")))
gmt_list2 <- lapply(gmt_list, function(x) list(x[1], x[2], x[3:length(x)]))

gmt_df <- as_tibble(do.call("rbind", gmt_list2)) %>%
  mutate(across(V1:V2, unlist)) %>%
  dplyr::rename(geneset=V1, link=V2, genes=V3) %>%
  filter(str_sub(geneset, 1, 4) == "GOBP") # filter to BP 

gmt_mat <- gmt_df %>%
  dplyr::select(-c(link)) %>%
  unnest(genes) %>%
  filter(genes %in% eb_distinct$EB_HGNC) %>%
  pivot_wider(names_from = genes, values_from = genes, values_fn = length, values_fill = 0) %>%
  column_to_rownames("geneset")

# subset to genesets with between 10 and 200 genes
valid_genesets <- (10 <= rowSums(gmt_mat)) & (rowSums(gmt_mat) < 200)
gmt_mat <- t(as.matrix(gmt_mat[which(valid_genesets),]))

is_novel <- !rownames(gmt_mat) %in% gtex_distinct$EB_HGNC
ora <- gseasusie::fit_ora(X=gmt_mat, y=is_novel)
ora$bhHypergeometric <- p.adjust(ora$pHypergeometric, method="BH")
ora$bhFishersExact <- p.adjust(ora$pFishersExact, method="BH")

gsea_results <- select(ora, c(geneSet, geneListSize, geneSetSize,
                              overlap, nGenes, propInList, propInSet,
                              oddsRatio, pFishersExact, bhFishersExact)) %>%
  write_tsv(gsea_table_loc)
