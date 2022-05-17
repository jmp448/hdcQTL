import os
import pandas as pd
from anndata import read_h5ad
import scdrs
import scanpy as sc

# read input files
project_loc = "/project2/gilad/jpopp/ebQTL"
adata_loc = snakemake.input[0]
gs_loc = snakemake.input[1]
gs_map_loc = snakemake.input[2]

# read output files
output_with_ctrl_loc = snakemake.output[0]
output_without_ctrl_loc = snakemake.output[1]

# read wildcards
geneset = snakemake.wildcards['geneset']
trait = snakemake.wildcards['trait']

# pull geneset and trait info from output file
results_loc = snakemake.output[0].rpartition("/")[0]
# geneset = results_loc.rpartition("/")[0].rpartition("/")[1]
intermediate_loc = os.path.join(results_loc, "intermediates/") 

# create folder if necessary at results/scDRS/{geneset}/{trait}/intermediate
if not os.path.exists(intermediate_loc):
    os.makedirs(intermediate_loc)

# load anndata object and gene set files
adata = read_h5ad(adata_loc)
df_gs = pd.read_csv(gs_loc, sep="\t")
trait_map = pd.read_csv(gs_map_loc, sep="\t")

# confirm expression data is log counts
adata.X = adata.layers['lognorm']

# compute stats
scdrs.preprocess(adata)

# create a dictionary mapping trait abbreviations to line entries in the gs file
trait_dict = dict(zip(trait_map['Abbreviation'], trait_map['Entry']))

# run scDRS
gene_list=df_gs[df_gs['TRAIT']==trait_dict[trait]]['GENESET'].values[0].split(",")
df_res = scdrs.score_cell(adata, gene_list, save_intermediate=intermediate_loc, return_ctrl_raw_score=True, return_ctrl_norm_score=True)

# save output
df_res.to_csv(output_with_ctrl_loc, sep="\t")
df_res.drop(df_res.filter(regex = 'ctrl').columns, axis = 1, inplace = True)
df_res.to_csv(output_without_ctrl_loc, sep="\t")
