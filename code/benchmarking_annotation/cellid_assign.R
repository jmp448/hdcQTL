library(Matrix)
library(Matrix.utils)
library(scran)
library(SingleCellExperiment)
library(tidyverse)
library(vroom)
library(caret)
library(CelliD)
library(corrplot)
library(reticulate)
use_condaenv("scvi-scanpy")
library(zellkonverter)

# Load cell id data
load("/project2/gilad/jpopp/ebQTL/data/fca/cellid_signatures_and_hvgs.Rdata")

# Load EB data
path = '/project2/gilad/jpopp/ebQTL/data/single_cell_objects/highpass/eb_pflog1ppfnorm.fca_hvgs.h5ad' 
eb <- zellkonverter::readH5AD(
  file=path,
  X_name="logcounts",
  use_hdf5=TRUE,
  uns=FALSE,
  varm=FALSE, obsm=FALSE,
  varp=FALSE, obsp=FALSE,
  verbose=T
)
# eb <- eb[intersect(hvg.fca_subsampled, rownames(eb)),] # filter to FCA's HVGs

eb <- RunMCA(eb, nmcs=50) # QUESTIONABLE
cellid_enrichments <- RunCellHGT(eb, pathways = new_signatures, dims = 1:50, n.features=100)

cellid_labels <- rownames(cellid_enrichments)[apply(cellid_enrichments, 2, which.max)]
cellid_labels <- ifelse(apply(cellid_enrichments, 2, max)>2, yes = cellid_labels, "unassigned") 

eb$celltype <- cellid_labels

saveRDS(eb, file="/project2/gilad/jpopp/ebQTL/data/single_cell_objects/highpass/eb_pflog1ppfnorm.fca_hvgs.mca.sce")
# writeH5AD(eb, file="/project2/gilad/jpopp/ebQTL/data/single_cell_objects/highpass/eb_pflog1ppfnorm.fca_hvgs.mca.h5ad", X_name = "logcounts", skip_assays = FALSE)

# write labels to a tsv file
as_tibble(cellid_labels, rownames="cell") %>%
  write_tsv("/project2/gilad/jpopp/ebQTL/data/fca/cellid_labels.tsv")

# write all assignment information to a tsv file
cellid_enrichments %>%
  as.matrix %>%
  t %>%
  as_tibble(rownames="cell") %>%
  write_tsv("/project2/gilad/jpopp/ebQTL/data/fca/cellid_enrichments.tsv")
