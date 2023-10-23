library(nullranges)
library(vroom)
library(tidyverse)
library(data.table)
set.seed(1234)

tests_loc <- snakemake@input[['tests_bed']]
hits_loc <- snakemake@input[['hits_bed']]
af_loc <- snakemake@input[['afs']]
gtf_loc <- snakemake@input[['gtf']]
bim_loc <- snakemake@input[['bim']]

background_bedfiles <- snakemake@output

# Load the full set of tests
all_tests <- vroom(tests_loc)

bim <- vroom(bim_loc, col_names=c("#CHR", "variant_id", "POS_CM", "POS_BP", "ALLELE_1", "ALLELE_2"),
      col_select=c("#CHR", "POS_BP", "variant_id")) %>%
  dplyr::rename(END=`POS_BP`) %>%
  mutate(START=as.integer(END)-1)

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
  mutate(ALT_FREQ=as.numeric(ALT_FREQ), REF_FREQ=as.numeric(REF_FREQ)) %>%
  mutate(MAF=map2_dbl(REF_FREQ, ALT_FREQ, min)) %>%
  inner_join(all_tests, by=c("CHROM"="#CHR", "POS"="END"), multiple="all") # allow for variants tested at multiple genes

all_tests <- all_tests %>%
  left_join(distinct(select(afs, c(EB_VARIANT_ID, MAF))), by="EB_VARIANT_ID") %>%
  drop_na() # drop any tests we didn't compute AF for

## Assign distance to TSS for each test
gencode <- vroom(gtf_loc) %>%
  rowwise() %>%
  mutate(tss=if_else(strand=="+", start, end))

all_tests <- all_tests %>%
  inner_join(select(gencode, c(ensg, tss)), by=c("EB_ENSG"="ensg")) %>%
  mutate(dist2tss=map2_dbl(tss, END, function(x,y){abs(x-y)}))

### Match the hits
hits <- vroom(hits_loc) %>%
  select(EB_VARIANT_ID, EB_HGNC) %>%
  filter(EB_VARIANT_ID != ".") %>%
  distinct() %>%
  inner_join(select(all_tests, c(EB_HGNC, EB_VARIANT_ID, MAF, dist2tss)), by=c("EB_HGNC", "EB_VARIANT_ID")) 

hits_dt <- hits %>%
  unite(qtl_test, EB_VARIANT_ID, EB_HGNC, sep="--") %>%
  dplyr::rename(qtl_MAF=MAF, qtl_TSS=dist2tss) %>%
  as.data.table

nonhits_dt <- select(all_tests, c(EB_HGNC, EB_VARIANT_ID, MAF, dist2tss)) %>%
  unite(bg_test, EB_VARIANT_ID, EB_HGNC, sep="--") %>%
  filter(!bg_test %in% hits_dt$qtl_test) %>%
  dplyr::rename(bg_MAF=MAF, bg_TSS=dist2tss) %>%
  as.data.table

# For each eQTL, list variants in the pool with MAF w/in 0.05 and gene distance w/in 5kb
setkeyv(hits_dt, c("qtl_MAF", "qtl_TSS"))
setkeyv(nonhits_dt, c("bg_MAF", "bg_TSS"))

# Define the thresholds
maf_window <- 0.05
tss_window <- 5000

hits_dt[, `:=`(MAF_lower = qtl_MAF - maf_window, MAF_upper = qtl_MAF + maf_window,
               TSS_lower = qtl_TSS - tss_window, TSS_upper = qtl_TSS + tss_window)]

# Map matched background tests to each QTL 
sample_sets <- nonhits_dt[hits_dt, on = .(bg_MAF >= MAF_lower, bg_MAF <= MAF_upper,
                                          bg_TSS >= TSS_lower, bg_TSS <= TSS_upper),
                          nomatch=0]

# Retain background tests' MAF, TSS distance, and now the QTLs it maps to
setkeyv(nonhits_dt, "bg_test")
simple_sample_sets <- nonhits_dt[sample_sets$bg_test]
simple_sample_sets$qtl_test <- sample_sets$qtl_test

# Convert to tibble for convenience moving forward
background_sets <- as_tibble(simple_sample_sets) %>%
  group_by(qtl_test) %>%
  sample_n(1000, replace=T) %>%
  group_by(qtl_test) %>%
  mutate(group_id = row_number()) %>%
  ungroup

background_bed_combined <- background_sets %>%
  select(bg_test, group_id) %>%
  separate(bg_test, into=c("variant_id", "hgnc"), sep="--") %>%
  left_join(bim, by=c("variant_id")) %>%
  left_join(distinct(select(gencode, c(hgnc, ensg))), by="hgnc") %>%
  dplyr::rename(EB_ENSG=ensg, EB_HGNC=hgnc, EB_VARIANT_ID=variant_id) %>%
  mutate(EB_TRAJECTORY="eb-cm") %>%
  relocate(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID, EB_TRAJECTORY, group_id)

# Save 1K background test sets
stopifnot(max(background_bed_combined$group_id) == length(background_bedfiles))
for(i in seq(length(background_bedfiles))) {
  write_tsv(background_bed_combined[background_bed_combined$group_id == i, -c(8)], background_bedfiles[[i]])
}
