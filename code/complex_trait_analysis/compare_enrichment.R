library(nullranges)
library(vroom)
library(tidyverse)
set.seed(1234)

# tests_list <- "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt"
# af_loc <- "data/genotypes/af_all.frq"
# mash_sighits_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash-signif_variant_gene_pairs.tsv"
# mash_sighits_overlap_loc <- "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/mash-signif_variant_gene_pairs.full_gtex_overlap.bed"
# eb_sighits_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.tsv"
# eb_sighits_overlap_loc <- "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.full_gtex_overlap.bed"
# gtf_loc <- "data/gencode/gencode.hg38.filtered.gtf"
# gtf_loc <- "/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
# gmt_file <- "data/gene_sets/c5.go.bp.v2022.1.Hs.symbols.gmt"
  
tests_list <- snakemake@input[['tests_list']]
af_loc <- snakemake@input[['afs']]
eb_sighits_loc <- snakemake@input[['eb_hits']]
eb_sighits_overlap_loc <- snakemake@input[['eb_gtex_overlap']]
gtf_loc <- snakemake@input[['gtf']]
gmt_file <- snakemake@input[['gmt']]

variant_candidates_info_loc <- snakemake@output[['candidate_info']]
variant_candidates_loc <- snakemake@output[['candidates']]
matchers_loc <- snakemake@output[['match_details']]

## List the genes with and without gtex overlap
gtex_overlap <- vroom(eb_sighits_overlap_loc, col_names=c("CHR", "START", "STOP", "EB_ENSG", "EB_HGNC", 
                                                          "RSID", "CELLTYPE", "CHR1", "START1", 
                                                          "STOP1", "GTEX_ENSG", "GTEX_REF", "GTEX_ALT"))
# If using mash-signif
eb_all_hits <- vroom(eb_sighits_loc) 
mash_all_hits <- vroom(mash_sighits_loc)
# If using signif
eb_all_hits <- vroom(eb_sighits_loc) %>%
  dplyr::rename(EB_HGNC=phenotype_id, EB_VARIANT_ID=variant_id)

overlap_egenes <- unique(gtex_overlap$EB_HGNC)
all_egenes <- unique(eb_all_hits$EB_HGNC)
novel_egenes <- setdiff(all_egenes, overlap_egenes)

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
eb_novel <- filter(eb_all_hits, EB_HGNC %in% novel_egenes) %>%
  dplyr::rename(phenotype_id_hgnc=EB_HGNC, variant_id_rsid=EB_VARIANT_ID) %>%
  inner_join(all_tests, by=c("phenotype_id_hgnc", "variant_id_rsid")) %>%
  select(c(variant_id, phenotype_id_ensg, MAF, dist2tss, variant_id_rsid)) %>% 
  distinct() %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  column_to_rownames("test")

eb_novel_eqtls <- as_tibble(eb_novel) %>%
  mutate(group="novel") 

## Matched background 1 - cis eQTLs with GTEx overlap
ciseqtl_matched <- filter(eb_all_hits, EB_HGNC %in% overlap_egenes) %>%
  dplyr::rename(phenotype_id_hgnc=EB_HGNC, variant_id_rsid=EB_VARIANT_ID) %>%
  inner_join(all_tests, by=c("phenotype_id_hgnc", "variant_id_rsid")) %>%
  select(c(variant_id, phenotype_id_ensg, MAF, dist2tss, variant_id_rsid)) %>%
  distinct() %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  column_to_rownames("test")

ciseqtl_matcher <- matchRanges(focal=eb_novel,
                               pool=ciseqtl_matched,
                               covar= ~ MAF + dist2tss,
                               method='stratified',
                               replace = FALSE)

ciseqtl_matched_eqtls <- as_tibble(ciseqtl_matched[ciseqtl_matcher@matchedIndex,]) %>%
  mutate(group="gtex_overlap")

## Matched background 3 - misc variants (non-eQTLs)
covariate_matched <- select(all_tests, c(variant_id, variant_id_rsid, phenotype_id_ensg, MAF, dist2tss)) %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  filter(!variant_id %in% c(eb_novel$variant_id, ciseqtl_matched$variant_id)) %>%
  column_to_rownames("test")

covariate_matcher <- matchRanges(focal=eb_novel,
                               pool=covariate_matched,
                               covar= ~ MAF + dist2tss,
                               method='stratified',
                               replace = FALSE)

covariate_matched_tests <- as_tibble(covariate_matched[covariate_matcher@matchedIndex,]) %>%
  mutate(group="covariate_matched")

all_followup_eqtls <- bind_rows(eb_novel_eqtls, 
                                ciseqtl_matched_eqtls,
                                covariate_matched_tests)

# Load schizophrenia summary stats
gwas_sumstats <- "data/gwas/zhang/sumstats/PASS_Schizophrenia_Pardinas2018.sumstats"
gwas_sumstats <- "data/gwas/zhang/sumstats/UKB_460K.body_HEIGHTz.sumstats"
gwas_sumstats <- "data/gwas/zhang/sumstats/UKB_460K.disease_CARDIOVASCULAR.sumstats"
gwas_sumstats <- "data/gwas/zhang/sumstats/UKB_460K.body_BMIz.sumstats"
gwas_sumstats <- "data/gwas/zhang/sumstats/UKB_460K.disease_HYPERTENSION_DIAGNOSED.sumstats"

scz <- vroom(gwas_sumstats)
p_norm <- 2*pnorm(q=abs(scz$Z), lower.tail=FALSE)
scz_shaky <- mutate(scz, p=2*pnorm(q=abs(scz$Z), lower.tail=FALSE))

scz_enrich <- inner_join(select(all_followup_eqtls, c(variant_id_rsid, group)), 
                         scz_shaky, by=c("variant_id_rsid"="SNP")) %>%
  arrange(group, p) %>% group_by(group) %>%
  add_tally() %>% mutate(group_index=row_number()) %>%
  mutate(nlog10_pval_unif=map2_dbl(group_index, n, function(a,b){-log10(a/b)})) %>%
  mutate(nlog10_pval=-log10(p)) %>%
  mutate(group=factor(group, levels=c("novel","gtex_overlap", "covariate_matched")))

ggplot(scz_enrich, aes(x = nlog10_pval_unif, y = nlog10_pval, color=group)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Expected (-log10 p-values)",
       y = "Observed (-log10 p-values)",
       title = "QQ Plot") +
  theme_minimal()

## We don't need to separate novel from GTEx overlap SNPs
eb_eqtls <- eb_all_hits %>%
  dplyr::rename(phenotype_id_hgnc=EB_HGNC, variant_id_rsid=EB_VARIANT_ID) %>%
  inner_join(all_tests, by=c("phenotype_id_hgnc", "variant_id_rsid")) %>%
  select(c(variant_id, phenotype_id_ensg, MAF, dist2tss, variant_id_rsid)) %>% 
  distinct() %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  column_to_rownames("test") %>%
  mutate(group="eb_eqtls")

mash_eqtls <- mash_all_hits %>%
  dplyr::rename(phenotype_id_hgnc=EB_HGNC, variant_id_rsid=EB_VARIANT_ID) %>%
  inner_join(all_tests, by=c("phenotype_id_hgnc", "variant_id_rsid")) %>%
  select(c(variant_id, phenotype_id_ensg, MAF, dist2tss, variant_id_rsid)) %>% 
  distinct() %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  column_to_rownames("test") %>%
  mutate(group="mash_eqtls")

covariate_matched <- select(all_tests, c(variant_id, variant_id_rsid, phenotype_id_ensg, MAF, dist2tss)) %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  filter(!variant_id %in% c(our_eqtls$variant_id)) %>%
  column_to_rownames("test")

covariate_matched_mash <- select(all_tests, c(variant_id, variant_id_rsid, phenotype_id_ensg, MAF, dist2tss)) %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  filter(!variant_id %in% c(mash_eqtls$variant_id)) %>%
  column_to_rownames("test")

covariate_matcher_ebeqtls <- matchRanges(focal=eb_eqtls,
                                 pool=covariate_matched,
                                 covar= ~ MAF + dist2tss,
                                 method='stratified',
                                 replace = FALSE)

covariate_matcher_masheqtls <- matchRanges(focal=mash_eqtls,
                                          pool=covariate_matched,
                                          covar= ~ MAF + dist2tss,
                                          method='stratified',
                                          replace = FALSE)

covariate_matched_toebeqtls <- as_tibble(covariate_matched[covariate_matcher_ebeqtls@matchedIndex,]) %>%
  mutate(group="covariate_matched")

covariate_matched_tomasheqtls <- as_tibble(covariate_matched[covariate_matcher_masheqtls@matchedIndex,]) %>%
  mutate(group="covariate_matched")

followup_eqtls_nosplit <- bind_rows(eb_eqtls, covariate_matched_toebeqtls)
followup_eqtls_nosplit_mash <- bind_rows(mash_eqtls, covariate_matched_tomasheqtls)

## Height
height_sumstats <- "data/gwas/zhang/sumstats/UKB_460K.body_HEIGHTz.sumstats"
height_gwas <- vroom(height_sumstats)
height_shaky <- mutate(height_gwas, p=2*pnorm(q=abs(height_gwas$Z), lower.tail=FALSE))

# Height - EB eQTLS
height_enrich <- inner_join(select(followup_eqtls_nosplit, c(variant_id_rsid, group)), 
                         height_shaky, by=c("variant_id_rsid"="SNP")) %>%
  arrange(group, p) %>% group_by(group) %>%
  add_tally() %>% mutate(group_index=row_number()) %>%
  mutate(nlog10_pval_unif=map2_dbl(group_index, n, function(a,b){-log10(a/b)})) %>%
  mutate(nlog10_pval=-log10(p)) %>%
  mutate(group=factor(group, levels=c("eb_eqtls","covariate_matched")))

ggplot(height_enrich, aes(x = nlog10_pval_unif, y = nlog10_pval, color=group)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Expected (-log10 p-values)",
       y = "Observed (-log10 p-values)",
       title = "QQ Plot") +
  theme_minimal() +
  ggtitle("Height, EB eQTLs")

# Height - MASH eQTLS
height_enrich_mash <- inner_join(select(followup_eqtls_nosplit_mash, c(variant_id_rsid, group)), 
                            height_shaky, by=c("variant_id_rsid"="SNP")) %>%
  arrange(group, p) %>% group_by(group) %>%
  add_tally() %>% mutate(group_index=row_number()) %>%
  mutate(nlog10_pval_unif=map2_dbl(group_index, n, function(a,b){-log10(a/b)})) %>%
  mutate(nlog10_pval=-log10(p)) %>%
  mutate(group=factor(group, levels=c("mash_eqtls","covariate_matched")))

ggplot(height_enrich_mash, aes(x = nlog10_pval_unif, y = nlog10_pval, color=group)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Expected (-log10 p-values)",
       y = "Observed (-log10 p-values)",
       title = "QQ Plot") +
  theme_minimal()

## SCZ
scz_sumstats <- "data/gwas/zhang/sumstats/PASS_Schizophrenia_Pardinas2018.sumstats"
scz_gwas <- vroom(scz_sumstats)
scz_shaky <- mutate(scz_gwas, p=2*pnorm(q=abs(scz_gwas$Z), lower.tail=FALSE))

### EB eQTLs
scz_enrich <- inner_join(select(followup_eqtls_nosplit, c(variant_id_rsid, group)), 
                            scz_shaky, by=c("variant_id_rsid"="SNP")) %>%
  arrange(group, p) %>% group_by(group) %>%
  add_tally() %>% mutate(group_index=row_number()) %>%
  mutate(nlog10_pval_unif=map2_dbl(group_index, n, function(a,b){-log10(a/b)})) %>%
  mutate(nlog10_pval=-log10(p)) %>%
  mutate(group=factor(group, levels=c("eb_eqtls","covariate_matched")))

ggplot(scz_enrich, aes(x = nlog10_pval_unif, y = nlog10_pval, color=group)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Expected (-log10 p-values)",
       y = "Observed (-log10 p-values)",
       title = "QQ Plot") +
  theme_minimal() +
  ggtitle("SCZ, EB eQTLs")

### MASH eQTLs
scz_enrich_mash <- inner_join(select(followup_eqtls_nosplit_mash, c(variant_id_rsid, group)), 
                         scz_shaky, by=c("variant_id_rsid"="SNP")) %>%
  arrange(group, p) %>% group_by(group) %>%
  add_tally() %>% mutate(group_index=row_number()) %>%
  mutate(nlog10_pval_unif=map2_dbl(group_index, n, function(a,b){-log10(a/b)})) %>%
  mutate(nlog10_pval=-log10(p)) %>%
  mutate(group=factor(group, levels=c("mash_eqtls","covariate_matched")))

ggplot(scz_enrich_mash, aes(x = nlog10_pval_unif, y = nlog10_pval, color=group)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Expected (-log10 p-values)",
       y = "Observed (-log10 p-values)",
       title = "QQ Plot") +
  theme_minimal()

