## Single-cell preprocessing and cell type annotation
### Alignment and demultiplexing
Alignment is at `/project2/gilad/kenneth/Pipelines/HumanCellranger/Snakefile_cellranger` with data in `kenneth/ebQTL/highpass/batch1/1001/*`
Katie's filtering and creation of the raw data is scattered between `/project2/gilad/katie/ebQTL/highpass_combinedFiles/102andPilot_MakeAdata.py`, `AddVireoAndFormColMetadata_AfterCellrangerAgr102andPilot.ipynb`, `102andPilot_filter.ipynb`, and `102andPilot_AddQC_ApplyFilter.ipynb`. 
These went into the creation of the raw_data file.
All of Kenneth's stuff he has, and Katie's stuff is copied into the `filtering` subfolder in `code` and `data`

### Classifier development (fetal cell atlas)
- Raw fetal cell atlas data can be downloaded from `https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/` (5,000 sampled cells from each cell type, cell annotations, and gene annotations)
- The workflow for generating fetal cell type gene signatures is at `rules/fca_annotation.py`. We explored three versions of the classifier: one with all 77 fetal cell atlas subtypes, one with a subset of 33 cell types curated to be more distinctive and found to offer greater classification accuracy on held-out test data, and one augmented set including the 33 fetal cell atlas cell types along with an IPSC signature from Plurinet
  - The relevant code is found in the `code/annotation` folder
  - Note that the 33 cell type classifier was used for the annotation of cells that was ultimately used for eQTL calling across cell types

### Cell type annotation
- The workflow for generating the fetal cell atlas classifier and classifying HDC cells is at `rules/fca_annotation.py`
- Code for this part of the analysis is in the folder `code/annotation`

### UMAP visualization and hierarchical clustering of cell types
- `rules/sc_preprocessing.py` contains the rules used to generate the UMAP embedding

## eQTL calling 
### Pseudobulk aggregation
- UMI counts for each donor and cell type are first aggregated in `analysis/annotation/assign_cellid.ipynb`
- Initial pseudobulk normalization and generation of a few QC metrics is done in `code/static_qtl_calling/pseudobulk_tmm-basic-qc.R`, followed by manual QC (based on inspection of PC plots) in `analysis/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/pseudobulk_qc.Rmd`. After removing outlier samples 
based on this manual step, we preprocess the pseudobulk data in `code/static_qtl_calling/pseudobulk_tmm-basic-agg.R`
- 

### eQTL calling in each cell type

### Mash multivariate eQTL analysis

## Interaction eQTL calling
### Trajectory isolation and pseudotime inference
### Dynamic eQTL calling
### Topic modeling
### Topic eQTL calling

## Software notes
- To run the `mash` analysis, we used `flashier`, which was not available on anaconda at the time
  - We therefore installed the package manually from https://github.com/willwerscheid/flashier 
  into the conda environment specified at `slurmy/r-mashr.yml` 



