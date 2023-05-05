library(tidyverse)

# Load the single cell data and filter to the cell types we're looking at in EBs
sc <- readRDS("../data/seurat.annotated.rds")
sc <- sc[,sc$type %in% c("IPSC", "MES", "CMES", "PROG", "CM")]
counts <- sc@assays[["RNA"]]@counts %>% t

# Save pseudobulk
type.ind <- paste(sc$individual, sc$type, sep="_")
counts.type <- counts %>% aggregate.Matrix(type.ind, fun="sum") %>% t %>% 
  as_tibble(rownames="gene") %>% arrange(gene) %>%
  write_tsv("data/static_qtl_calling/elorbany_cmstages/pseudobulk_tmm/elorbany_cmstages.pseudobulk_tmm.tsv")

# Get sample summary
sample_summary <- tibble(ind_type=paste(sc$individual, sc$type, sep="_"), libsize=sc$nCount_RNA) %>%
  mutate(ncells=1) %>%
  group_by(ind_type) %>% summarize(total_counts=sum(libsize), n_cells_unfiltered=sum(ncells)) %>%
  mutate(individual=str_extract(ind_type, "[^_]+"), type=str_extract(ind_type, "[^_]+$"), n_cells_filtered=n_cells_unfiltered) %>%
  mutate(dropped=total_counts < 1e5) %>%
  relocate(ind_type, individual, type, n_cells_unfiltered, dropped, total_counts, n_cells_filtered) %>%
  write_tsv("data/static_qtl_calling/elorbany_cmstages/pseudobulk_tmm/sample_summary.tsv")
