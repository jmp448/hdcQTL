library(Matrix)
library(Matrix.utils)
library(scran)
library(SingleCellExperiment)
library(tidyverse)
library(vroom)
library(caret)
library(CelliD)
library(corrplot)

# Load subsampled fetal cell atlas data
fca_counts <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/counts.sampled.rds")
fca_cells <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/cell_metadata.rds")
fca_genes <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/gene_metadata.rds") %>%
  mutate_all(as.character)

# Subset cells to those in the subsampled data
fca_cells <- fca_cells %>%
  filter(sample %in% colnames(fca_counts))

# Subset genes to protein-coding genes without duplicated HGNC entries
pull_gene_type <- function(attr) {
  str_split(attr, "\"")[[1]][6]
}
pull_gene_name <- function(attr) {
  str_split(attr, "\"")[[1]][8]
}
pull_gene_id <- function(attr) {
  str_split(attr, "\"")[[1]][2]
}

gencode <- vroom("/project2/gilad/kenneth/References/human/cellranger/cellranger4.0/refdata-gex-GRCh38-2020-A/genes/genes.gtf", 
                 col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=5) %>%
  filter(seqname %in% paste0("chr", seq(1, 22))) %>%
  filter(feature == "gene") %>%
  mutate(type=map_chr(attribute, pull_gene_type)) %>%
  mutate(ensg=map_chr(attribute, pull_gene_id)) %>%
  mutate(hgnc=map_chr(attribute, pull_gene_name)) %>%
  filter(type=="protein_coding")

fca_genes <- mutate(fca_genes, gene_id_short=str_extract(gene_id, "[^.]+"))
fca_pc_genes <- filter(fca_genes, 
                       (gene_short_name %in% gencode$hgnc) & (gene_id_short %in% gencode$ensg))
hgnc_duplicates <- fca_pc_genes %>%
  group_by(gene_short_name) %>%
  count() %>%
  filter(n > 1) %>%
  pull(gene_short_name)
fca_pc_genes <- fca_pc_genes %>%
  filter(!gene_short_name %in% hgnc_duplicates)

# Make the counts matrix match the metadata
fca_counts <- fca_counts[fca_pc_genes$gene_id,fca_cells$sample]
rownames(fca_counts) <- fca_pc_genes$gene_short_name

# Create SingleCellExperiment object
fca <- SingleCellExperiment(list(counts=fca_counts),
                            colData=fca_cells)
saveRDS(fca, "/project2/gilad/jpopp/ebQTL/data/fca/counts.sampled.sce")

# Normalize expression data
fca <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/counts.sampled.sce")

# Subsample each of the major cell types to 5000 or the total number of cells
fca$celltype <- str_extract(fca$Organ_cell_lineage, "[^-]+$")
celltype_sampler <- as_tibble(colData(fca)) %>%
  select(sample, celltype) %>%
  group_by(celltype) %>%
  add_count(celltype)
keepers_rare <- celltype_sampler %>%
  filter(n <= 5000) %>%
  pull(sample)
keepers_sampled <- celltype_sampler %>%
  filter(n > 5000) %>%
  slice_sample(n=5000) %>%
  pull(sample)
fca <- fca[,c(keepers_rare, keepers_sampled)]

libs <- librarySizeFactors(fca)
sizeFactors(fca) <- libs
fca <- logNormCounts(fca)

fca.dec <- modelGeneVar(fca)
hvg.fca <- getTopHVGs(fca.dec, fdr.threshold=0.1)

# Dimensionality reduction
### PCA
fca.pca <- prcomp(t(logcounts(fca)[hvg.fca,]), rank=50)
reducedDims(fca) <- list(PCA=fca.pca$x)

### MCA
fca <- RunMCA(fca, nmcs=50, features=hvg.fca)
saveRDS(fca, "/project2/gilad/jpopp/ebQTL/data/fca/counts.sampled.mca_pca_fromhvgs.sce")

# Check the performance of a Cell ID classifier on the dataset it came from (FCA)
fca_signatures <- GetGroupGeneSet(fca, group.by="celltype", n.features=100)

## Look at which gene sets are tightly overlapping


## Test the classification scheme performance on the original data
cellid_enrichments <- RunCellHGT(fca, pathways = fca_signatures, features=hvg.fca, dims = 1:50, n.features=100)

cellid_labels <- rownames(cellid_enrichments)[apply(cellid_enrichments, 2, which.max)]
cellid_labels <- ifelse(apply(cellid_enrichments, 2, max)>2, yes = fca_labels, "unassigned") #6 unassigned cells here

### View confusion matrix
major_types <- c(unique(colData(fca)$celltype), "unassigned")
confusion <- confusionMatrix(data=factor(cellid_labels, levels=major_types), reference = factor(colData(fca)$celltype, levels=major_types))

corrplot(confusion$table, is.corr=F, tl.pos='lt', tl.cex=0.5)


# Modify which cell types to include
fca <- readRDS("/project2/gilad/jpopp/ebQTL/data/fca/counts.sampled.sce")
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
  
# Update the dataset to contain even samples from these groups
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

libs <- librarySizeFactors(fca_subsampled)
sizeFactors(fca_subsampled) <- libs
fca_subsampled <- logNormCounts(fca_subsampled)

fca_subsampled.dec <- modelGeneVar(fca_subsampled)
hvg.fca_subsampled <- getTopHVGs(fca_subsampled.dec, fdr.threshold=0.1)

tibble(gene=hvg.fca_subsampled) %>%
  write_tsv("/project2/gilad/jpopp/ebQTL/data/fca/fca_subsampled_hvg.tsv")

# Dimensionality reduction
### PCA
fca_subsampled.pca <- prcomp(t(logcounts(fca_subsampled)[hvg.fca_subsampled,]), rank=50)
reducedDims(fca_subsampled) <- list(PCA=fca_subsampled.pca$x)

### MCA
fca_subsampled <- RunMCA(fca_subsampled, nmcs=50, features=hvg.fca_subsampled)

# Check the performance of a Cell ID classifier on the dataset it came from (FCA)
new_signatures <- GetGroupGeneSet(fca_subsampled, group.by="celltype", n.features=100)

new_cellid_enrichments <- RunCellHGT(fca_subsampled, pathways = new_signatures, features=hvg.fca_subsampled, dims = 1:50, n.features=100)

new_cellid_labels <- rownames(new_cellid_enrichments)[apply(new_cellid_enrichments, 2, which.max)]
new_cellid_labels <- ifelse(apply(new_cellid_enrichments, 2, max)>2, yes = new_cellid_labels, "unassigned") #3 unassigned cells here

### View confusion matrix
major_types_remaining <- c(unique(colData(fca_subsampled)$celltype), "unassigned")
confusion_temp <- confusionMatrix(data=factor(new_cellid_labels, levels=major_types_remaining), reference = factor(colData(fca_subsampled)$celltype, levels=major_types_remaining))

corrplot(confusion_temp$table, is.corr=F, tl.pos='lt', tl.cex=0.75, method='color')

# Cell type tree
mca.embedding.subsampled <- reducedDim(fca_subsampled, "MCA")
new_celltype_centroids <- aggregate.Matrix(mca.embedding.subsampled, groupings=fca_subsampled$celltype, fun="mean")

new_dists <- dist(new_celltype_centroids)
hc <- hclust(new_dists)
plot(hc, cex=0.8)

# Visualize MC space
emb <- as_tibble(reducedDim(fca_subsampled, "MCA"), rownames="sample") %>%
  left_join(as_tibble(colData(fca_subsampled)), by="sample") 
centroids <- as_tibble(as.matrix(new_celltype_centroids), rownames="celltype")

ggplot(emb, aes(x=MCA_1, y=MCA_2, color=celltype)) +
  geom_point() +
  geom_label(data=centroids, aes(label=celltype))

# Visualize distance in terms of hamming distance between cell types
celltypes <- unique(names(new_signatures))
ct_dists <- matrix(nrow=length(celltypes), ncol=length(celltypes))
rownames(ct_dists) <- celltypes
colnames(ct_dists) <- celltypes
for (ct1 in seq_len(length(celltypes))) {
  for (ct2 in seq_len(length(celltypes))) {
    ct_dists[ct1, ct2] = 100 - length(intersect(new_signatures[[ct1]], new_signatures[[ct2]]))
  }
}
hc <- hclust(as.dist(ct_dists))
plot(hc, cex=0.8)

celltypes <- unique(names(fca_signatures))
ct_dists <- matrix(nrow=length(celltypes), ncol=length(celltypes))
rownames(ct_dists) <- celltypes
colnames(ct_dists) <- celltypes
for (ct1 in seq_len(length(celltypes))) {
  for (ct2 in seq_len(length(celltypes))) {
    ct_dists[ct1, ct2] = 100 - length(intersect(fca_signatures[[ct1]], fca_signatures[[ct2]]))
  }
}
hc <- hclust(as.dist(ct_dists))
plot(hc, cex=0.8)

# Save cell type signatures, highly variable genes, and our filtered object
saveRDS(fca_subsampled, "/project2/gilad/jpopp/ebQTL/data/fca/counts.subsampled.mca_pca_fromhvgs.sce")
save(new_signatures, hvg.fca_subsampled, file="/project2/gilad/jpopp/ebQTL/data/fca/cellid_signatures_and_hvgs.Rdata")
  