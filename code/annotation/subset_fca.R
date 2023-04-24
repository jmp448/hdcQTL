library(Matrix)
library(Matrix.utils)
library(scran)
library(SingleCellExperiment)
library(tidyverse)
set.seed(42)

# Read snakemake inputs
fca_counts_loc <- snakemake@input[['counts']]
fca_cell_metadata_loc <- snakemake@input[['cell_metadata']]
fca_gene_metadata_loc <- snakemake@input[['gene_metadata']]
pcg_loc <- snakemake@input[['pc_genes']]
subsampled_loc <- snakemake@output[['fca_subsampled']]

# Load subsampled fetal cell atlas data
fca_counts <- readRDS(fca_counts_loc)
fca_cells <- readRDS(fca_cell_metadata_loc)
fca_genes <- readRDS(fca_gene_metadata_loc) %>%
  mutate_all(as.character)

pc_genes <- read_tsv(pcg_loc)

# Subset cells to those in the subsampled data
fca_cells <- fca_cells %>%
  filter(sample %in% colnames(fca_counts))

# Subset genes to protein-coding genes without duplicated HGNC entries
rownames(fca_counts) <- str_extract(rownames(fca_counts), "[^.]+")
fca_pc_genes <- pc_genes %>%
  filter(ensg %in% rownames(fca_counts)) %>%
  arrange(ensg)

# Make the counts matrix match the metadata
fca_counts <- fca_counts[fca_pc_genes$ensg, fca_cells$sample] # subset genes & cells
rownames(fca_counts) <- fca_pc_genes$hgnc # switch genes to hgnc

# Create SingleCellExperiment object
fca <- SingleCellExperiment(list(counts=fca_counts),
                            colData=fca_cells)

# Revise cell type labels
fca$celltype <- str_extract(fca$Organ_cell_lineage, "[^-]+$")

### Remove any placenta cell types, or cell types that overlapped another cell type or were uncharacterized
omitted_types <- c("AFP_ALB positive cells",
                   "CCL19_CCL21 positive cells",
                   "CLC_IL5RA positive cells",
                   "CSH1_CSH2 positive cells",
                   "ELF3_AGBL2 positive cells",
                   "MUC13_DMBT1 positive cells",
                   "PDE11A_FAM19A2 positive cells",
                   "PDE1C_ACSM3 positive cells",
                   "SATB2_LRRC7 positive cells",
                   "SKOR2_NPSR1 positive cells",
                   "SLC24A4_PEX5L positive cells",
                   "SLC26A4_PAEP positive cells")

celltype_sampler <- as_tibble(colData(fca)) %>%
  select(sample, Organ, celltype) %>%
  mutate(celltype_revised=celltype) %>%
  mutate(celltype_revised=if_else(Organ=="Placenta", true=NA_character_, false=celltype_revised)) %>%
  mutate(celltype_revised=if_else(celltype %in% omitted_types, true=NA_character_, false=celltype_revised))

### BLOOD CELLS
### Remove antigen presenting cells as there aren't enough of them to get a distinct signature
### Merge lymphoid cells and thymocytes
### Remove myeloid cells and hematopoietic stem/ progenitor cells since they're ancestors of other cell types
### Keep microglia on their own
ancestor_types <- c("Myeloid cells", "Hematopoietic stem cells")

celltype_sampler <- celltype_sampler %>%
  mutate(celltype_revised=if_else(celltype == "Antigen presenting cells", true=NA_character_, false=celltype_revised)) %>%
  mutate(celltype_revised=if_else(celltype %in% ancestor_types, true=NA_character_, false=celltype_revised)) %>%
  mutate(celltype_revised=if_else(celltype == "Thymocytes", true="Lymphoid cells", false=celltype_revised))

### EYE
### Remove lens fibre cells since there are only ~200, corneal and conjunctival epithelial cells since ~100
### Merge retinal progenitors & Muller glia, amacrine cells, bipolar cells, ganglion cells, retinal pigment cells, horizontal cells, photoreceptor cells under retinal cells
retinal_types <- c("Retinal progenitors and Muller glia", "Amacrine cells",
                   "Bipolar cells", "Ganglion cells", "Retinal pigment cells",
                   "Horizontal cells", "Photoreceptor cells")

celltype_sampler <- celltype_sampler %>%
  mutate(celltype_revised=if_else(celltype == "Lens fibre cells", true=NA_character_, false=celltype_revised)) %>%
  mutate(celltype_revised=if_else(celltype == "Corneal and conjunctival epithelial cells", true=NA_character_, false=celltype_revised)) %>%
  mutate(celltype_revised=if_else(celltype %in% retinal_types, true="Retinal cells", false=celltype_revised))

### CNS, PNS, Neuroendocrine CONSOLIDATION
cns_neuronal_types <- c("Limbic system neurons", "Excitatory neurons",
                        "Inhibitory neurons", "Purkinje neurons",
                        "Granule neurons", "Inhibitory interneurons",
                        "Unipolar brush cells")
cns_glial_types <- c("Microglia", "Astrocytes", "Oligodendrocytes")
pns_neuronal_types <- c("Visceral neurons", "ENS neurons")
pns_glial_types <- c("ENS glia", "Satellite cells", "Schwann cells")
neuroendocrine_types <- c("Chromaffin cells", "Islet endocrine cells",
                          "Neuroendocrine cells", "Sympathoblasts")

celltype_sampler <- celltype_sampler %>%
  mutate(celltype_revised=if_else(celltype %in% cns_neuronal_types, true="CNS neurons", false=celltype_revised)) %>%
  mutate(celltype_revised=if_else(celltype %in% cns_glial_types, true="CNS glia", false=celltype_revised)) %>%
  mutate(celltype_revised=if_else(celltype %in% pns_neuronal_types, true="PNS neurons", false=celltype_revised)) %>%
  mutate(celltype_revised=if_else(celltype %in% pns_glial_types, true="PNS glia", false=celltype_revised)) %>%
  mutate(celltype_revised=if_else(celltype %in% neuroendocrine_types, true="Neuroendocrine cells", false=celltype_revised))
  
# Finally, remove thymic epithelial cells, it's the last somewhat rare cell population
celltype_sampler <- celltype_sampler %>%
  mutate(celltype_revised=if_else(celltype == "Thymic epithelial cells", true=NA_character_, false=celltype_revised))
  
# Update the dataset to contain even samples from these cell types, and a max of 5k cells per major cell type (regardless of organ)
celltype_subsetter <- celltype_sampler %>%
  drop_na() %>%
  group_by(celltype_revised) %>%
  add_count(celltype_revised)
keepers_rare <- celltype_subsetter %>%
  filter(n <= 5000) %>%
  pull(sample)
keepers_sampled <- celltype_subsetter %>%
  filter(n > 5000) %>%
  slice_sample(n=5000) %>%
  pull(sample)
celltype_subsetter <- celltype_subsetter %>%
  filter(sample %in% c(keepers_rare, keepers_sampled)) %>%
  select(-celltype) %>%
  rename(celltype=celltype_revised)

fca_subsampled <- fca[,celltype_subsetter$sample]
fca_subsampled$celltype <- celltype_subsetter$celltype

saveRDS(fca_subsampled, subsampled_loc)