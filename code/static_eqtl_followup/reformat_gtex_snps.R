library(tidyverse)
library(vroom)

# snp_loc <- "/project2/gilad/jpopp/GTEx_Analysis_v8_eQTL/Adipose_Subcutaneous.v8.signif_variant_gene_pairs.txt"

options(future.globals.maxSize=200*1000^3)
snp_loc <- snakemake@input[[1]]
bed_loc <- snakemake@output[[1]]

snps <- vroom(snp_loc, col_select=c("variant_id", "gene_id")) 
snps <- mutate(snps, GTEX_GENE=str_extract(gene_id, "[^.]+"), .keep="unused")
snps <- separate(snps, variant_id, into=c("#CHR", "END", "GTEX_REF", "GTEX_ALT", "build"), sep="_")
snps <- mutate(snps, START=as.numeric(END)-1)
snps <- select(snps, c(`#CHR`, START, END, GTEX_GENE, GTEX_REF, GTEX_ALT)) 
write_tsv(snps, bed_loc)
