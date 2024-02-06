library(tidyverse)

fca_signatures_loc <- snakemake@input[['fca_signatures']]
tidy_signatures_loc <- snakemake@output[['fca_signatures_tsv']]

# Create a table of gene set signatures, save to tsv
fca_signature_genesets <- readRDS(fca_signatures_loc)

fca_genesets_tsv <- as_tibble(fca_signature_genesets) %>% 
  as.data.frame() %>% t %>% as_tibble(rownames="celltype") %>%
  pivot_longer(!celltype, names_to="temp", values_to="gene") %>%
  select(-c(temp)) %>%
  write_tsv(tidy_signatures_loc)