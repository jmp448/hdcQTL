library(tidyverse)
library(vroom)

# ldblock_loc <- "results/static_eqtl_followup/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/yri_tagging.list.geno.ld"
# tophits_loc <- "results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.pos"
# gtf_loc <- "/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

ldblock_loc <- snakemake@input[['ld_block']]
tophits_loc <- snakemake@input[['top_hits']]
gtf_loc <- snakemake@input[['gtf_loc']]
bim_file <- snakemake@input[['bim_file']]

ldblock_bed_loc <- snakemake@output[[1]]

yri_ld <- vroom(ldblock_loc)
eqtls <- vroom(tophits_loc)
bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS", "variant_id"))

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
  filter(type=="protein_coding")

# 4 HGNC symbols are duplicated here - these have overlapping loci, and will be removed from our analysis
# TBCE, MATR3, HSPA14, GGT1
hgnc.dup <- filter(gencode, hgnc %in% names(table(gencode$hgnc)[table(gencode$hgnc)>1]))
gencode <- gencode %>% filter(!hgnc %in% hgnc.dup$hgnc)
hgnc_ensg_dict <- select(gencode, c(hgnc, ensg))

eqtl_block <- left_join(yri_ld, eqtls, by=c("CHR1"="#CHR", "POS1"="POS"), multiple="all") %>%
  left_join(bim, by=c("CHR1"="#CHR", "POS1"="POS")) %>%
  pivot_longer(cols = c(CHR1, POS1, CHR2, POS2),
               names_to = c(".value", "set"),
               names_pattern = "(\\w+)(\\d)") %>%
  select(CHR, POS, GENE, variant_id) %>%
  dplyr::rename(`#CHR`=CHR, END=POS, EB_VARIANT_ID=variant_id) %>%
  mutate(START=END-1, .before=END) %>%
  left_join(hgnc_ensg_dict, by=c("GENE"="hgnc")) %>%
  dplyr::rename(EB_ENSG=ensg, EB_HGNC=GENE) %>%
  relocate(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID) %>%
  write_tsv(ldblock_bed_loc)
