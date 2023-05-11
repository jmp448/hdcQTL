#! /usr/bin/env Rscript

# Install package 'anndata' to read anndata (h5ad)
# Install package 'sceasy' to convert h5ad to seurat object (.rds)
# Required modules: R, gcc - 'module load R gcc'
# Installation: devtools::install_github("cellgeni/sceasy"); ref: https://github.com/cellgeni/sceasy 

# Load a few packages.
library(sceasy)

input <- snakemake@input[['anndata']]
outfile <- snakemake@output[['seurat']]

sceasy::convertFormat(input, from="anndata", to="seurat", outFile=outfile)


#  shell:
#         """
#         Rscript /project2/gilad/mli/eb_project/snakemake/2.h5ad2seurat/h5ad2seurat.R -i {input.input_file} --from {params.from_} --to {params.to_} -o {output.output_file}
#         """
