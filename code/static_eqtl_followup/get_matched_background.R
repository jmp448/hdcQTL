library(nullranges)
library(vroom)
library(tidyverse)
set.seed(1234)

tests_loc <- snakemake@input[['tests_bed']]
hits_loc <- snakemake@input[['hits_bed']]
af_loc <- snakemake@input[['afs']]
gtf_loc <- snakemake@input[['gtf']]

matchable_hits_output <- snakemake@output[['hits_bed']]
matched_background_output <- snakemake@output[['background_bed']]
matcher_loc <- snakemake@output[['match_details']]

# Load the full set of tests
all_tests <- vroom(tests_loc)

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

hits_df <- hits %>%
  unite(test, EB_VARIANT_ID, EB_HGNC, sep="-", remove=F) %>%
  column_to_rownames("test")

nonhits_df <- select(all_tests, c(EB_HGNC, EB_VARIANT_ID, MAF, dist2tss)) %>%
  unite(test, EB_VARIANT_ID, EB_HGNC, sep="-", remove=F) %>%
  filter(!test %in% rownames(hits_df)) %>%
  column_to_rownames("test")

matcher <- matchRanges(focal=hits_df, 
                  pool=nonhits_df, 
                  covar= ~ MAF + dist2tss,
                  method='stratified',
                  replace = FALSE)
saveRDS(matcher, matcher_loc)

matched_hits <- as_tibble(nonhits_df[matcher@matchedIndex,]) %>%
  left_join(all_tests, by=c("EB_VARIANT_ID", "EB_HGNC")) %>%
  select(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID) %>%
  write_tsv(matched_background_output)

# Also save a new subset of the significant hits which we were able to include in the focal set
matchable_hits <- as_tibble(hits_df) %>%
  left_join(all_tests, by=c("EB_VARIANT_ID", "EB_HGNC")) %>%
  select(`#CHR`, START, END, EB_ENSG, EB_HGNC, EB_VARIANT_ID) %>%
  write_tsv(matchable_hits_output)
