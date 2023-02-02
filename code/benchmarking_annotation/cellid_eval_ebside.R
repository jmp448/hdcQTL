library(tidyverse)
library(Matrix.utils)
library(fgsea)
library(susieR)
library(pathways)
library(vroom)
library(gseasusie)

### Load & tidy data
# Read in the input files from Snakefile
eb_loc <- snakemake@input[["fca"]]
eb_loc <- "/project2/gilad/jpopp/ebQTL/data/single_cell_objects/highpass/eb_pflog1ppfnorm.fca_hvgs.mca.h5ad"
celltype_signatures <- snakemake@input[["signatures"]]
table_prefix <- snakemake@params[['table_prefix']]

# Load cell and gene embeddings
eb <- readRDS("/project2/gilad/jpopp/ebQTL/data/single_cell_objects/highpass/eb_pflog1ppfnorm.fca_hvgs.mca.sce")
eb$celltype[eb$celltype=="unassigned"] <- NA
# Get a list of cell types
cell.types <- unique(eb$celltype)
cell.types <- cell.types[!is.na(cell.types)]

# Pull MCA embeddings of cells and genes
cell_emb <- as_tibble(reducedDim(eb, "MCA"), rownames="cell") %>%
  left_join(as_tibble(colData(eb), rownames="cell"), by="cell") 
gene_emb <- attributes(cell_emb$MCA_1)$genesCoordinates

# Get cell type centroids
celltype_centroids <- as.matrix(aggregate.Matrix(reducedDim(eb, "MCA"), groupings=eb$celltype, fun="mean"))

### Distribution of cells around the cell type centroid
# For each cell type, measure the distance from each cell to the centroid
list_cell_dists <- function(t) {
  cell_emb_t <- cell_emb %>%
    filter(celltype==t) %>%
    select(c(cell, starts_with("MCA"))) %>%
    column_to_rownames("cell") %>%
    as.matrix()
  centroid_t <- t(celltype_centroids[t,])
  diffs_t <- sweep(cell_emb_t, 2, centroid_t)
  dists_t <- apply(diffs_t, 1, norm, type="2") %>%
    as_tibble() %>%
    add_column(type=t)
  return(dists_t)
}

cell_dists <- map_dfr(cell.types, list_cell_dists)
write_tsv(cell_dists, "/project2/gilad/jpopp/ebQTL/data/fca/eb_centroid_cell_dists.tsv")

ggplot(filter(cell_dists, value<= 10), aes(x=value, fill=type)) + 
  geom_density() +
  facet_wrap(vars(type))

### Distribution of celltype marker genes around the cell type centroid
# For each cell type, measure the distance from each gene to the centroid
list_gene_dists <- function(t) {
  gene_emb_t <- gene_emb[intersect(rownames(gene_emb), new_signatures[[t]]),]
  centroid_t <- t(celltype_centroids[t,])
  diffs_t <- sweep(gene_emb_t, 2, centroid_t)
  dists_t <- apply(diffs_t, 1, norm, type="2") %>%
    as_tibble(rownames="gene") %>%
    add_column(type=t)
  return(dists_t)
}

gene_dists <- map_dfr(cell.types, list_gene_dists)
write_tsv(gene_dists, "/project2/gilad/jpopp/ebQTL/data/fca/eb_centroid_gene_dists.tsv")

ggplot(gene_dists, aes(x=value, fill=type)) + 
  geom_density() +
  facet_wrap(vars(type))

# Focus on some of the more surprising cell types 
# For lymphoid and megakaryocyte cell types, find 100 closest genes
find_closest <- function(t, num.neighbors=100) {
  centroid_t <- t(celltype_centroids[t,])
  diffs_t <- sweep(gene_emb, 2, centroid_t)
  markers_t <- apply(diffs_t, 1, norm, type="2") %>%
    as_tibble(rownames="gene") %>%
    arrange(value) %>% slice_head(n=num.neighbors) %>%
    pull(gene)
  return(markers_t)
}

lymphoid_closest_genes <- lapply(c("Lymphoid cells"), find_closest)
lymphoid_drivers <- intersect(lymphoid_closest_genes[[1]], new_signatures[['Lymphoid cells']])

megakaryocyte_closest_genes <- lapply(c("Megakaryocytes"), find_closest)
megakaryocyte_drivers <- intersect(megakaryocyte_closest_genes[[1]], new_signatures[['Megakaryocytes']])

# For each of the epithelial cell types, find their 100 closest genes
epithelial_celltypes <- c("Acinar cells", 
                          "Metanephric cells", 
                          #"Ciliated epithelial cells",
                          #"Intestinal epithelial cells", 
                          #"Bronchiolar and alveolar epithelial cells",
                          #"Ureteric bud cells",
                          "Ductal cells")


epithelial_closest_genes <- lapply(epithelial_celltypes, find_closest)
overlap_epithelial_genes <- Reduce(intersect, epithelial_closest_genes)
# Scroll down for GSEA on this group of genes



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

## Pathway enrichment analysis - lymphoid annotation driver genes
y_lymphoid <- bg.hgnc %in% lymphoid_drivers

logistic.fit <- gseasusie::fit_logistic_susie_veb_boost(X_pathways, y_lymphoid, L=10)
ora <- gseasusie::fit_ora(X_pathways, y_lymphoid)

if (!is.null(logistic.fit$sets$cs)) {
  cs <- pivot_longer(as_tibble(logistic.fit$sets$cs), everything(), names_to="component", values_to="index")
} else {
  cs <- tibble(component=character(), index=numeric())
}

credible_sets <- as_tibble(logistic.fit$pip, rownames="geneSet") %>%
  dplyr::rename(pip=value) %>%
  rowid_to_column("index") %>%
  left_join(cs, by="index") %>%
  left_join(ora, by="geneSet")

gsea_pathways_lymphoid <- credible_sets %>%
  left_join(select(genesets_pathways_info, c(name, id, description_brief)), by=c("geneSet"="id")) %>%
  relocate(name, pip, component, pFishersExact)

### Megakaryocyte driver genes
## Pathway enrichment analysis
y_megakaryocyte <- bg.hgnc %in% megakaryocyte_drivers

logistic.fit.mega.pathway <- gseasusie::fit_logistic_susie_veb_boost(X_pathways, y_megakaryocyte, L=10)
ora.mega.pathway <- gseasusie::fit_ora(X_pathways, y_megakaryocyte)

if (!is.null(logistic.fit.mega.pathway$sets$cs)) {
  cs_mega.mega.pathway <- pivot_longer(as_tibble(logistic.fit.mega.pathway$sets$cs), everything(), names_to="component", values_to="index")
} else {
  cs <- tibble(component=character(), index=numeric())
}

credible_sets_mega_pathway <- as_tibble(logistic.fit.mega.pathway$pip, rownames="geneSet") %>%
  dplyr::rename(pip=value) %>%
  rowid_to_column("index") %>%
  left_join(cs, by="index") %>%
  left_join(ora, by="geneSet")

gsea_pathways_megakaryocyte <- credible_sets_mega_pathway %>%
  left_join(select(genesets_pathways_info, c(name, id, description_brief)), by=c("geneSet"="id")) %>%
  relocate(name, pip, component, pFishersExact)

## Celltype enrichment analysis
logistic.fit.mega.celltype <- gseasusie::fit_logistic_susie_veb_boost(X_celltype, y_megakaryocyte, L=10)
ora.mega.celltype <- gseasusie::fit_ora(X_celltype, y_megakaryocyte)

if (!is.null(logistic.fit.mega.pathway$sets$cs)) {
  cs_mega.mega.pathway <- pivot_longer(as_tibble(logistic.fit.mega.pathway$sets$cs), everything(), names_to="component", values_to="index")
} else {
  cs <- tibble(component=character(), index=numeric())
}

credible_sets_mega_pathway <- as_tibble(logistic.fit.mega.pathway$pip, rownames="geneSet") %>%
  dplyr::rename(pip=value) %>%
  rowid_to_column("index") %>%
  left_join(cs, by="index") %>%
  left_join(ora, by="geneSet")

gsea_pathways_megakaryocyte <- credible_sets_mega_pathway %>%
  left_join(select(genesets_pathways_info, c(name, id, description_brief)), by=c("geneSet"="id")) %>%
  relocate(name, pip, component, pFishersExact)


## Pathway enrichment analysis - epithelial overlap genes
y_epi <- bg.hgnc %in% pull(filter(acinar_drivers, value < 5), "gene")
  
logistic.fit <- gseasusie::fit_logistic_susie_veb_boost(X_celltype, y_epi, L=10)
ora <- gseasusie::fit_ora(X_celltype, y_epi)

cs <- pivot_longer(as_tibble(logistic.fit$sets$cs), everything(), names_to="component", values_to="index")

credible_sets <- as_tibble(logistic.fit$pip, rownames="geneSet") %>%
  rename(pip=value) %>%
  rowid_to_column("index") %>%
  left_join(cs, by="index") %>%
  left_join(ora, by="geneSet") 

gsea_celltype_epi <- credible_sets %>%
  left_join(select(genesets_celltype_info, c(name, id, description_brief)), by=c("geneSet"="id")) %>%
  relocate(name, pip, component, pFishersExact)

susie_gsea_pathways <- gsea_pathways %>%
  drop_na(component)

write_tsv(gsea_pathways, "/project2/gilad/jpopp/ebQTL/data/fca/cellid_genesets_pathways_gsea.tsv")
