library(tidyverse)
library(Matrix.utils)
library(fgsea)
library(susieR)
library(pathways)
library(vroom)
library(gseasusie)

### Load & tidy data
# Read in the input files from Snakefile
fca_loc <- snakemake@input[["fca"]]
celltype_signatures <- snakemake@input[["signatures"]]
table_prefix <- snakemake@params[['table_prefix']]

# Load cell and gene embeddings
fca <- readRDS(fca_loc)

# Get a list of cell types
cell.types <- unique(fca$celltype)

# Load cell type signatures
load(celltype_signatures)

# Pull MCA embeddings of cells and genes
cell_emb <- as_tibble(reducedDim(fca, "MCA"), rownames="sample") %>%
  left_join(as_tibble(colData(fca)), by="sample") 
gene_emb <- attributes(cell_emb$MCA_1)$genesCoordinates

# Get cell type centroids
celltype_centroids <- as.matrix(aggregate.Matrix(reducedDim(fca, "MCA"), groupings=fca$celltype, fun="mean"))

### Distribution of cells around the cell type centroid
# For each cell type, measure the distance from each cell to the centroid
list_cell_dists <- function(t) {
  cell_emb_t <- cell_emb %>%
    filter(celltype==t) %>%
    select(c(sample, starts_with("MCA"))) %>%
    column_to_rownames("sample") %>%
    as.matrix()
  centroid_t <- t(celltype_centroids[t,])
  diffs_t <- sweep(cell_emb_t, 2, centroid_t)
  dists_t <- apply(diffs_t, 1, norm, type="2") %>%
    as_tibble() %>%
    add_column(type=t)
  return(dists_t)
}

cell_dists <- map_dfr(cell.types, list_cell_dists)
write_tsv(cell_dists, "/project2/gilad/jpopp/ebQTL/data/fca/cellid_centroid_cell_dists.tsv")

ggplot(filter(cell_dists, value<= 10), aes(x=value, fill=type)) + 
  geom_density() +
  facet_wrap(vars(type))

### Distribution of genes around the cell type centroid
# For each cell type, measure the distance from each gene to the centroid
list_gene_dists <- function(t) {
  gene_emb_t <- gene_emb[new_signatures[[t]],]
  centroid_t <- t(celltype_centroids[t,])
  diffs_t <- sweep(gene_emb_t, 2, centroid_t)
  dists_t <- apply(diffs_t, 1, norm, type="2") %>%
    as_tibble() %>%
    add_column(type=t)
  return(dists_t)
}

gene_dists <- map_dfr(cell.types, list_gene_dists)
write_tsv(gene_dists, "/project2/gilad/jpopp/ebQTL/data/fca/cellid_centroid_gene_dists.tsv")

ggplot(gene_dists, aes(x=value, fill=type)) + 
  geom_density() +
  facet_wrap(vars(type))

### GSEA for each marker gene set
# For each cell type, perform GSEA on the marker gene set
# Translate the FCA's HGNC names to ENSG
gtf_loc <- "/project2/gilad/jpopp/ebQTL/data/gencode/gencode.hg38.filtered.gtf"
pull_gene_type <- function(attr) {
  str_split(attr, "\"")[[1]][11]
}
pull_gene_name <- function(attr) {
  str_split(attr, "\"")[[1]][15]
}
pull_gene_id <- function(attr) {
  str_split(attr, "\"")[[1]][3]
}
gencode <- vroom(gtf_loc, col_names=c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"), skip=5) %>%
  filter(seqname %in% paste0("chr", seq(1, 22))) %>%
  filter(feature == "gene") %>%
  mutate(type=map_chr(attribute, pull_gene_type)) %>%
  mutate(hgnc=map_chr(attribute, pull_gene_name)) %>%
  mutate(ensg=map_chr(attribute, pull_gene_id)) %>%
  filter(type=="protein_coding")

## Cell type enrichment analysis
genesets_celltype <- gene_sets_human$gene_set_info %>%
  filter(database == "MSigDB-C8") %>%
  pull(id)
genesets_celltype_info <- gene_sets_human$gene_set_info %>%
  filter(database=="MSigDB-C8")

# Background gene list (in both HGNC and ENSG format)
bg.ensg <- gencode %>%
  filter(hgnc %in% hvg.fca_subsampled) %>%
  pull(ensg) %>%
  intersect(rownames(gene_sets_human$gene_sets)) %>%
  sort
bg.hgnc <- gencode %>%
  filter(ensg %in% bg.ensg) %>%
  arrange(ensg) %>%
  pull(hgnc)

X_celltype <- gene_sets_human$gene_sets[bg.ensg, genesets_celltype]
X_celltype <- X_celltype[,colSums(X_celltype) >= 10]

logistic_gsea <- function(t) {
  y_t <- bg.hgnc %in% new_signatures[[t]]
  
  logistic.fit <- gseasusie::fit_logistic_susie_veb_boost(X_celltype, y_t, L=10)
  ora <- gseasusie::fit_ora(X_celltype, y_t)
  
  if (!is.null(logistic.fit$sets$cs)) {
    cs <- pivot_longer(as_tibble(logistic.fit$sets$cs), everything(), names_to="component", values_to="index")
  } else {
    cs <- tibble(component=character(), index=numeric())
  }
  
  credible_sets <- as_tibble(logistic.fit$pip, rownames="geneSet") %>%
    rename(pip=value) %>%
    rowid_to_column("index") %>%
    left_join(cs, by="index") %>%
    left_join(ora, by="geneSet") %>%
    mutate(fca_type=!!t)
}

gsea_celltypes <- map_dfr(cell.types, logistic_gsea) %>%
  left_join(select(genesets_celltype_info, c(name, id, description_brief)), by=c("geneSet"="id")) %>%
  relocate(name, fca_type, pip, component, pFishersExact)

susie_gsea_celltypes <- gsea_celltypes %>%
  drop_na(component)

write_tsv(gsea_celltypes, "/project2/gilad/jpopp/ebQTL/data/fca/cellid_genesets_celltype_gsea.tsv")

## Pathway enrichment analysis
genesets_pathways <- gene_sets_human$gene_set_info %>%
  filter(database == "BioSystems-kegg") %>%
  pull(id)
genesets_pathways_info <- gene_sets_human$gene_set_info %>%
  filter(database=="BioSystems-kegg")

X_pathways <- gene_sets_human$gene_sets[bg.ensg, genesets_pathways]
X_pathways <- X_pathways[,colSums(X_pathways) >= 5]

logistic_gsea_pathways <- function(t) {
  y_t <- bg.hgnc %in% new_signatures[[t]]
  
  logistic.fit <- gseasusie::fit_logistic_susie_veb_boost(X_pathways, y_t, L=10)
  ora <- gseasusie::fit_ora(X_pathways, y_t)
  
  if (!is.null(logistic.fit$sets$cs)) {
    cs <- pivot_longer(as_tibble(logistic.fit$sets$cs), everything(), names_to="component", values_to="index")
  } else {
    cs <- tibble(component=character(), index=numeric())
  }
  
  credible_sets <- as_tibble(logistic.fit$pip, rownames="geneSet") %>%
    rename(pip=value) %>%
    rowid_to_column("index") %>%
    left_join(cs, by="index") %>%
    left_join(ora, by="geneSet") %>%
    mutate(fca_type=!!t)
}

gsea_pathways <- map_dfr(cell.types, logistic_gsea_pathways) %>%
  left_join(select(genesets_pathways_info, c(name, id, description_brief)), by=c("geneSet"="id")) %>%
  relocate(name, fca_type, pip, component, pFishersExact)

susie_gsea_pathways <- gsea_pathways %>%
  drop_na(component)

write_tsv(gsea_pathways, "/project2/gilad/jpopp/ebQTL/data/fca/cellid_genesets_pathways_gsea.tsv")
