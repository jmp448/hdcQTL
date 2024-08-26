library(tidyverse)

scdrs_celltype_loc <- snakemake@input
output_loc <- snakemake@output[[1]]

celltypes_combined <- read_tsv(snakemake@input, id="trait") %>%
mutate(trait=basename(dirname(trait))) %>%
dplyr::rename(celltype=group) %>%
write_tsv(output_loc)

