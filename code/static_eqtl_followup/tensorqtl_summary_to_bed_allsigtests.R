library(tidyverse)
library(vroom)

eqtl_file="results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/signif_variant_gene_pairs.tsv"
bim_file="data/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/Retinal-cells/genotypes_filtered_plink.bim"
gtf_loc="/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
#bed_file="results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.bed"

eqtl_file <- snakemake@input[['qtl_summary']]
bim_file <- snakemake@input[['bim_file']]
gtf_loc <- snakemake@input[['gtf_loc']]

bed_file <- snakemake@output[['bedfile']]

## Get a mapping from rsID to chromosome positions
bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS_BP", "variant_id")) %>%
  dplyr::rename(END=`POS_BP`) %>%
  mutate(START=as.integer(END)-1) # bed files are zero-indexed

# Get a mapping from HGNC to ENSG
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

# Make BED file
bed <- vroom(eqtl_file) %>%
  dplyr::rename(GENE=phenotype_id) %>%
  select(c(variant_id, GENE, celltype)) %>%
  left_join(bim, by="variant_id") %>%
  left_join(hgnc_ensg_dict, by=c("GENE"="hgnc")) %>%
  group_by(`#CHR`, START, END, ensg, GENE, variant_id) %>%
  summarize(celltype=paste(celltype, collapse=",")) %>%
  dplyr::rename(EB_ENSG=ensg, EB_HGNC=GENE, EB_VARIANT_ID=variant_id, EB_CELLTYPE=celltype) %>%
  relocate(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID, EB_CELLTYPE) %>%
  drop_na() %>%
  write_tsv(bed_file)
  