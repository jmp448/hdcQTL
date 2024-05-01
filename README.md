# Code for eQTL analyses in heterogeneous differentiating cultures 
## Single-cell preprocessing 
- Code for filtering is available in `code/filtering/`
- The output of this sequence is the raw data file available for download from GEO

### Classifier development (fetal cell atlas)
- Raw fetal cell atlas data can be downloaded from `https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/` (5,000 sampled cells from each cell type, cell annotations, and gene annotations)
- The workflow for generating fetal cell type gene signatures is at `rules/fca_annotation.py`. We explored three versions of the classifier: one with all 77 fetal cell atlas subtypes, one with a subset of 33 cell types curated to be more distinctive and found to offer greater classification accuracy on held-out test data, and one augmented set including the 33 fetal cell atlas cell types along with an IPSC signature from Plurinet
  - The relevant code is found in the `code/annotation` folder
  - Note that the 33 cell type classifier was used for the annotation of cells that was ultimately used for eQTL calling across cell types

### Cell type annotation
- The workflow for generating the fetal cell atlas classifier and classifying HDC cells is at `rules/fca_annotation.py`
- Code for this part of the analysis is in the folder `code/annotation`

### UMAP visualization and hierarchical clustering of cell types
- `rules/sc_preprocessing.py` contains the rules used to generate the UMAP embedding, with associated code in `code/sc_preprocessing/`

## Static eQTL calling 
### Pseudobulk aggregation
- UMI counts for each donor and cell type are first aggregated in `analysis/annotation/assign_cellid.ipynb`
- Initial pseudobulk normalization and generation of a few QC metrics is done in `code/static_qtl_calling/pseudobulk_tmm-basic-qc.R`, followed by manual QC (based on inspection of PC plots) in `analysis/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/pseudobulk_qc.Rmd`. After removing outlier samples based on this manual step, we preprocess the pseudobulk data in `code/static_qtl_calling/pseudobulk_tmm-basic-agg.R`

### eQTL calling in each cell type
- After data wrangling (spread over several rules in `rules/static_qtl_calling.py`), we call eQTLs using tensorQTL in `code/static_qtl_calling/static_qtl_calling/tensorqtl_permutations.py`
- We control the FDR across all genes in `code/static_qtl_calling/static_qtl_calling/tensorqtl_fdr.py`

### Mash multivariate eQTL analysis
- Snakemake rules for the multivariate (cross-celltype) mash analysis are in `rules/mash_qtl_calling.py`, with associated code in `code/mash_qtl_calling/`
## Interaction eQTL calling
### Trajectory isolation and pseudotime inference
- Snakemake rules for trajectory isolation and pseudotime estimation are in `rules/trajectory_inference.py`, with associated code in `code/trajectory_inference/`
### Dynamic eQTL calling
- Snakemake rules for pseudobulk aggregation (by pseudotime binning) and dynamic eQTL calling are in `rules/dynamic_qtl_calling.py`, with associated code in `code/dynamic_qtl_calling/`
### Topic modeling
- Snakemake rules for pseudocell aggregation, topic modeling, and the topic DE analysis are in `rules/fast_topics.py`, with associated code in `code/fast_topics/`
### Topic eQTL calling
- Snakemake rules for topic eQTL calling with CellRegMap, including multiple testing correction and *post hoc* estimation of effect sizes, are in `rules/cellregmap_eqtl_calling.py`, with associated code in `code/cellregmap_eqtl_calling/` 
## Complex trait analysis
- Snakemake rules for the comparison of cell type eQTLs to schizophrenia GWAS loci are in `rules/scz_analysis.py` with code in `
## Software notes
- To run the `mash` analysis, we used `flashier`, which was not available on anaconda at the time
  - We therefore installed the package manually from https://github.com/willwerscheid/flashier 
  into the conda environment specified at `slurmy/r-mashr.yml` 
