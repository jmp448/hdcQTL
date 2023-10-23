library(tidyverse)
library(vroom)
library(qvalue)
library(MatrixGenerics)

genotype_loc <- snakemake@input[['genotypes']]
qtls_loc <- snakemake@input[['qtls']]
pseudotime_loc <- snakemake@input[['pseudotime']]
bim_file <- snakemake@input[['bim_file']]
gtf_loc <- snakemake@input[['gtf_loc']]
traj <- snakemake@wildcards[['trajectory']]

early_bedfile <- snakemake@output[['early']]
late_bedfile <- snakemake@output[['late']]
combined_bedfile <- snakemake@output[['combined']]

classify_dynamics <- function(b_g, b_gi, tmin, tmax, gmin, gmax, threshold=1) {
  early_effect = (gmin*b_g + gmin*tmin*b_gi) - (gmax*b_g + gmax*tmin*b_gi)
  late_effect = (gmin*b_g + gmin*tmax*b_gi) - (gmax*b_g + gmax*tmax*b_gi)
  
  if (sign(early_effect) == sign(late_effect)) {
    if (abs(early_effect) >= abs(late_effect)) {
      class = "early"
    } else {
      class = "late"
    }
    # QTL switching signs
  } else {
    if ((abs(early_effect) >= abs(late_effect)) & (abs(late_effect) < threshold)) {
      class = "early"
    } else if ((abs(early_effect) < abs(late_effect)) & (abs(early_effect) < threshold)) {
      class = "late"
      # Switch magnitude is greater than or equal to threshold - is this switching behavior significant enough
    } else if ((abs(early_effect) >= threshold) & (abs(late_effect) >= threshold)) {
      class = "switch"
    } 
    else {
      class = "unassigned"
    }
  }
  class
}

genotypes <- vroom(genotype_loc) %>%
  column_to_rownames("snp")
gmins <- rowMins(as.matrix(genotypes))
gmaxs <- rowMaxs(as.matrix(genotypes))
granges <- tibble(variant_id=rownames(genotypes), gmin=gmins, gmax=gmaxs) 

emt <- vroom(qtls_loc)
pseudotime <- vroom(pseudotime_loc)

qtls <- filter(emt, pval_adj_bh <= 0.1) %>%
  left_join(granges, by="variant_id") %>%
  mutate(tmin=min(pseudotime$pseudotime), tmax=max(pseudotime$pseudotime)) %>%
  mutate(class=pmap_chr(.l=list(b_g, b_gi, tmin, tmax, gmin, gmax), .f=classify_dynamics))

# Convert to BED format
bim <- vroom(bim_file, col_names=c("#CHR", "variant_id", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2"),
             col_select=c("#CHR", "POS_BP", "variant_id")) %>%
  dplyr::rename(END=`POS_BP`) %>%
  mutate(START=as.integer(END)-1) # bed files are zero-indexed
gencode <- vroom(gtf_loc) %>% 
  select(c(hgnc, ensg))

qtls_bed <- qtls %>%
  filter(variant_id != ".") %>%
  dplyr::rename(GENE=phenotype_id) %>%
  select(c(variant_id, GENE, class)) %>%
  distinct() %>%
  inner_join(bim, by="variant_id") %>%
  inner_join(gencode, by=c("GENE"="hgnc")) %>%
  dplyr::rename(EB_ENSG=ensg, EB_HGNC=GENE, EB_VARIANT_ID=variant_id) %>%
  relocate(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID, class) %>%
  mutate(EB_TRAJECTORY=traj)

early_bed <- filter(qtls_bed, class %in% c("early")) %>%
  select(-c(class)) %>%
  write_tsv(early_bedfile)

late_bed <- filter(qtls_bed, class %in% c("late")) %>%
  select(-c(class)) %>%
  write_tsv(late_bedfile)