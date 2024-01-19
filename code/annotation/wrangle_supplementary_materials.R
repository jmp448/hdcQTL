library(tidyverse)

# Create a table of gene set signatures, save to tsv
fca_signature_genesets <- readRDS("data/fca/fca_signatures.rds")

fca_genesets_tsv <- as_tibble(fca_signature_genesets) %>% 
  as.data.frame() %>% t %>% as_tibble(rownames="celltype") %>%
  pivot_longer(!celltype, names_to="temp", values_to="gene") %>%
  select(-c(temp)) %>%
  write_tsv("data/fca/fca_signatures.tsv")