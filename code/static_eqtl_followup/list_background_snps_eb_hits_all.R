library(nullranges)
library(vroom)
library(tidyverse)
set.seed(1234)

# tests_list <- "results/static_eqtl_followup/eb_cellid/pseudobulk_tmm/basic/8pcs/eb_gtex_harmonized_tests.txt"
# af_loc <- "data/genotypes/af_all.frq"
# gtf_loc <- "/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
# # mash_hits_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash_sighits_lfsr_0.05.tsv"
# all_hits_loc <- "results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.tsv"
  
tests_list <- snakemake@input[['test_list']]
af_loc <- snakemake@input[['afs']]
gtf_loc <- snakemake@input[['gtf']]
all_hits_loc <- snakemake@input[['all_hits']]

matchable_hits_output <- snakemake@output[['hits_bed']]
matched_background_output <- snakemake@output[['background_bed']]
matcher_loc <- snakemake@output[['match_details']]

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
pull_gene_type <- function(attr) {
  str_split(attr, "\"")[[1]][6]
}

pull_gene_name <- function(attr) {
  str_split(attr, "\"")[[1]][8]
}

pull_gene_ensg <- function(attr) {
  str_split(attr, "\"")[[1]][2]
}

gencode <- vroom(gtf_loc, col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=5) %>%
  filter(seqname %in% paste0("chr", seq(1, 22))) %>%
  filter(feature == "gene") %>%
  mutate(type=map_chr(attribute, pull_gene_type)) %>%
  mutate(hgnc=map_chr(attribute, pull_gene_name)) %>%
  mutate(ensg=map_chr(attribute, pull_gene_ensg)) %>%
  filter(type=="protein_coding") %>%
  rowwise() %>%
  mutate(tss=if_else(strand=="+", start, end))

all_tests <- all_tests %>%
  mutate(variant_pos=as.numeric(str_extract(variant_id, "(?<=_)[^_]+(?=_.*$)"))) %>%
  inner_join(select(gencode, c(ensg, tss)), by=c("phenotype_id_ensg"="ensg")) %>%
  mutate(dist2tss=map2_dbl(tss, variant_pos, function(x,y){abs(x-y)}))

### Match
# mash_hits <- vroom(mash_hits_loc) %>%
#   inner_join(all_tests, by=c("EB_HGNC"="phenotype_id_hgnc", "EB_RSID"="variant_id_rsid")) %>%
#   select(c(variant_id, phenotype_id_ensg, MAF, dist2tss)) %>%
#   unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
#   column_to_rownames("test")
# 
# non_mash_hits <- select(all_tests, c(variant_id, phenotype_id_ensg, MAF, dist2tss)) %>%
#   unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
#   filter(!test %in% rownames(mash_hits)) %>%
#   column_to_rownames("test")

all_hits <- vroom(all_hits_loc) %>%
  dplyr::rename(phenotype_id_hgnc=phenotype_id, variant_id_rsid=variant_id) %>%
  inner_join(all_tests, by=c("phenotype_id_hgnc", "variant_id_rsid")) %>%
  select(c(variant_id, phenotype_id_ensg, MAF, dist2tss)) %>% 
  distinct() %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  column_to_rownames("test")

non_hits <- select(all_tests, c(variant_id, phenotype_id_ensg, MAF, dist2tss)) %>%
  unite(test, variant_id, phenotype_id_ensg, sep="-", remove=F) %>%
  filter(!test %in% rownames(all_hits)) %>%
  column_to_rownames("test")

matcher <- matchRanges(focal=all_hits, 
                  pool=non_hits, 
                  covar= ~ MAF + dist2tss,
                  method='stratified',
                  replace = FALSE)
saveRDS(matcher, matcher_loc)

matched_hits <- as_tibble(non_hits[matcher@matchedIndex,]) %>%
  left_join(all_tests, by=c("variant_id", "phenotype_id_ensg")) %>%
  separate(variant_id, into=c("CHROM", "POS", "REF", "ALT", "BUILD"), sep="_") %>%
  mutate(START=as.numeric(POS)-1) %>%
  dplyr::rename(`#CHR`=CHROM, END=POS, EB_ENSG=phenotype_id_ensg, EB_HGNC=phenotype_id_hgnc, 
         EB_VARIANT_ID=variant_id_rsid) %>%
  select(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID) %>%
  mutate(EB_CELLTYPE=NA_character_) %>%
  write_tsv(matched_background_output)

# Also save a new subset of the significant hits which we were able to include in the focal set
matchable_hits <- as_tibble(all_hits) %>%
  left_join(all_tests, by=c("variant_id", "phenotype_id_ensg")) %>%
  separate(variant_id, into=c("CHROM", "POS", "REF", "ALT", "BUILD"), sep="_") %>%
  mutate(START=as.numeric(POS)-1) %>%
  dplyr::rename(`#CHR`=CHROM, END=POS, EB_ENSG=phenotype_id_ensg, EB_HGNC=phenotype_id_hgnc, 
                EB_VARIANT_ID=variant_id_rsid) %>%
  select(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID) %>%
  mutate(EB_CELLTYPE=NA_character_) %>%
  write_tsv(matchable_hits_output)
