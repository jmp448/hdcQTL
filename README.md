## Single-cell preprocessing and cell type annotation
### Alignment and demultiplexing
Alignment is at `/project2/gilad/kenneth/Pipelines/HumanCellranger/Snakefile_cellranger`
with data in `kenneth/ebQTL/highpass/batch1/1001/*`
Katie's filtering and creation of the raw data is scattered between `/project2/gilad/katie/ebQTL/highpass_combinedFiles/102andPilot_MakeAdata.py`,
`AddVireoAndFormColMetadata_AfterCellrangerAgr102andPilot.ipynb`, `102andPilot_filter.ipynb`, and `102andPilot_AddQC_ApplyFilter.ipynb`. 
These went into the creation of the raw_data file.
All of Kenneth's stuff he has, and Katie's stuff is copied into the `filtering` subfolder in `code` and `data`

### Classifier development (fetal cell atlas)
- Raw fetal cell atlas data can be downloaded from `https://atlas.brotmanbaty.org/bbi/human-gene-expression-during-development/` (5,000 sampled cells from each cell type, cell annotations, and gene annotations)
- 

### Cell type annotation

### UMAP visualization and hierarchical clustering of cell types

## eQTL calling 
### eQTL calling in each cell type

### Mash multivariate eQTL analysis

## Interaction eQTL calling
### Trajectory isolation and pseudotime inference
### Dynamic eQTL calling
### Topic modeling
### Topic eQTL calling




