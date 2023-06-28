library(nullranges)
library(vroom)
library(tidyverse)
set.seed(1234)

# tests_list <- "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt"
# af_loc <- "data/genotypes/af_all.frq"
# eb_sighits_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.tsv"
# eb_sighits_overlap_loc <- "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.full_gtex_overlap.bed"
# gtf_loc <- "data/gencode/gencode.hg38.filtered.gtf"
# gtf_loc <- "/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
# gmt_file <- "data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
# 
# go_geneset <- "GOBP_TISSUE_DEVELOPMENT"
  
tests_list <- snakemake@input[['tests_list']]
af_loc <- snakemake@input[['afs']]
eb_sighits_loc <- snakemake@input[['eb_hits']]
eb_sighits_overlap_loc <- snakemake@input[['eb_gtex_overlap']]
gtf_loc <- snakemake@input[['gtf']]
gmt_file <- snakemake@input[['gmt']]

go_geneset <- snakemake@wildcards[['gs']]

variant_candidates_info_loc <- snakemake@output[['candidate_info']]
variant_candidates_loc <- snakemake@output[['candidates']]
matchers_loc <- snakemake@output[['match_details']]

## List the genes with and without gtex overlap
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
novel_egenes_gs <- intersect(gs_genes, setdiff(all_egenes, overlap_egenes))
gtex_egenes_gs <- intersect(gs_genes, overlap_egenes)

## Load tests
all_tests <- vroom(tests_list)

## Assign MAF to each test
### Note - we're computing MAF on the full panel, so if you look at the distribution of MAFs, 
#          it doesn't cutoff at 0.1 bc MAF will vary based on which donors were included
#          in that specific celltype analysis. This also means some variants we tested aren't
#          represented here, bc their MAF on the full panel falls below 0.05 even though it's over 
#          0.1 among the donors tested in a specific cell type
afs <- vroom(af_loc, col_names=c("CHROM", "POS", "N_ALLELES", "N_CHR", "REF_FREQ", "ALT_FREQ"),
             skip=1) %>%
  separate(REF_FREQ, into=c("REF", "REF_FREQ"), sep=":") %>%
  separate(ALT_FREQ, into=c("ALT", "ALT_FREQ"), sep=":") %>%
  mutate(variant_id=paste(CHROM, POS, REF, ALT, "b38", sep="_")) %>%
  filter(variant_id %in% all_tests$variant_id) %>%
  mutate(ALT_FREQ=as.numeric(ALT_FREQ), REF_FREQ=as.numeric(REF_FREQ)) %>%
  mutate(MAF=map2_dbl(REF_FREQ, ALT_FREQ, min))

all_tests <- all_tests %>%
  inner_join(select(afs, c(variant_id, MAF)), by="variant_id")

## Assign distance to TSS for each test
#### Again, for simplicity we're just going to pull gene locations for a 
# get a mapping from HGNC to ENSG
# pull_gene_type <- function(attr) {
#   # str_split(attr, "\"")[[1]][11] for the online version
#   str_split(attr, "\"")[[1]][6] # for the kenneth version
# }
# 
# pull_gene_name <- function(attr) {
#   # str_split(attr, "\"")[[1]][15] for the online version
#   str_split(attr, "\"")[[1]][8] # for the kenneth version
# }
# 
# pull_gene_ensg <- function(attr) {
#   # str_split(attr, "\"")[[1]][3] for the online version
#   str_split(attr, "\"")[[1]][2] # for the kenneth version
# }

gencode <- vroom(gtf_loc, col_select=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "type", "hgnc", "ensg")) %>%
  mutate(tss=if_else(strand=="+", start, end))

all_tests <- all_tests %>%
  mutate(variant_pos=as.numeric(str_extract(variant_id, "(?<=_)[^_]+(?=_.*$)"))) %>%
  inner_join(select(gencode, c(ensg, tss)), by=c("phenotype_id_ensg"="ensg")) %>%
  mutate(dist2tss=map2_dbl(tss, variant_pos, function(x,y){abs(x-y)}))

## List the top QTLs associated with the 'novel' eGenes (that also were tested in GTEx)
pathway_novel <- filter(eb_all_hits, phenotype_id %in% novel_egenes_gs) %>%
  dplyr::rename(phenotype_id_hgnc=phenotype_id, variant_id_rsid=variant_id) %>%
  inner_join(all_tests, by=c("phenotype_id_hgnc", "variant_id_rsid")) %>%
  select(c(variant_id, phenotype_id_ensg, MAF, dist2tss)) %>% 
  distinct() %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  column_to_rownames("test")

pathway_novel_candidates <- as_tibble(pathway_novel) %>%
  mutate(group="pathway_novel") 

# Save just this set of variants 
trans_qtl_variant_candidates <- pathway_novel_candidates %>%
  write_tsv(variant_candidates_info_loc)

trans_qtl_variant_candidates_ids <- trans_qtl_variant_candidates %>%
  select(variant_id) %>%
  write_tsv(variant_candidates_loc, col_names=F)

### BEFORE UNCOMMENTING - note that matching needs to change to accommodate the use of entire LD blocks (not just matching variants any more)
# ## Matched background 1 - tests also in the gene set, but did have GTEx overlap
# pathway_matched <- filter(eb_all_hits, phenotype_id %in% gtex_egenes_gs) %>%
#   dplyr::rename(phenotype_id_hgnc=phenotype_id, variant_id_rsid=variant_id) %>%
#   inner_join(all_tests, by=c("phenotype_id_hgnc", "variant_id_rsid")) %>%
#   arrange(pval_nominal) %>%
#   group_by(phenotype_id_hgnc) %>%
#   slice_head(n=1) %>%
#   ungroup() %>%
#   select(c(variant_id, phenotype_id_ensg, MAF, dist2tss)) %>% 
#   distinct() %>%
#   unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
#   filter(!variant_id %in% pathway_novel$variant_id) %>%
#   column_to_rownames("test")
# 
# pathway_matcher <- matchRanges(focal=pathway_novel, 
#                        pool=pathway_matched, 
#                        covar= ~ MAF + dist2tss,
#                        method='stratified',
#                        replace = FALSE)
# 
# pathway_matched_candidates <- as_tibble(pathway_matched[pathway_matcher@matchedIndex,]) %>%
#   mutate(group="pathway_matched") 
# 
# ## Matched background 2 - misc cis eQTLs with GTEx overlap
# ciseqtl_matched <- filter(eb_all_hits, phenotype_id %in% overlap_egenes) %>%
#   dplyr::rename(phenotype_id_hgnc=phenotype_id, variant_id_rsid=variant_id) %>%
#   inner_join(all_tests, by=c("phenotype_id_hgnc", "variant_id_rsid")) %>%
#   arrange(pval_nominal) %>%
#   group_by(phenotype_id_hgnc) %>%
#   slice_head(n=1) %>%
#   ungroup() %>%
#   select(c(variant_id, phenotype_id_ensg, MAF, dist2tss)) %>% 
#   distinct() %>%
#   unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
#   filter(!variant_id %in% c(pathway_novel$variant_id, pathway_matched$variant_id)) %>%
#   column_to_rownames("test")
# 
# ciseqtl_matcher <- matchRanges(focal=pathway_novel, 
#                                pool=ciseqtl_matched, 
#                                covar= ~ MAF + dist2tss,
#                                method='stratified',
#                                replace = FALSE)
# 
# ciseqtl_matched_candidates <- as_tibble(ciseqtl_matched[ciseqtl_matcher@matchedIndex,]) %>%
#   mutate(group="ciseqtl_matched") 
# 
# ## Matched background 3 - misc variants (non-eQTLs)
# covariate_matched <- select(all_tests, c(variant_id, phenotype_id_ensg, MAF, dist2tss)) %>%
#   unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
#   filter(!variant_id %in% c(pathway_novel$variant_id, pathway_matched$variant_id, ciseqtl_matched$variant_id)) %>%
#   column_to_rownames("test")
# 
# covariate_matcher <- matchRanges(focal=pathway_novel, 
#                                pool=covariate_matched, 
#                                covar= ~ MAF + dist2tss,
#                                method='stratified',
#                                replace = FALSE)
# 
# covariate_matched_candidates <- as_tibble(covariate_matched[covariate_matcher@matchedIndex,]) %>%
#   mutate(group="covariate_matched")
# 
# # Save all four sets of variants
# trans_qtl_variant_candidates <- bind_rows(pathway_novel_candidates,
#                                           pathway_matched_candidates,
#                                           ciseqtl_matched_candidates,
#                                           covariate_matched_candidates) %>%
#   write_tsv(variant_candidates_info_loc)
# 
# trans_qtl_variant_candidates_ids <- trans_qtl_variant_candidates %>%
#   select(variant_id) %>%
#   write_tsv(variant_candidates_loc, col_names=F)
# 
# save(pathway_matcher, ciseqtl_matcher, covariate_matcher, file=matchers_loc)
#   
#   
