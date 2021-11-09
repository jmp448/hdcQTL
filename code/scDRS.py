import os
import pandas as pd
from anndata import read_h5ad
import scdrs
import scanpy as sc

DATA_PATH = "/project2/gilad/jpopp/ebQTL/data/"
adata = read_h5ad(os.path.join(DATA_PATH, "single_cell_objects/Lowpass.3seqbatches.merged.scvi_processed_and_scaled.h5ad"))
df_gs = pd.read_csv(os.path.join(DATA_PATH, "scDRS/gs_file/magma_10kb_1000.74_traits.gs"), sep="\t")

adata.X = adata.layers['lognorm']

scdrs.method.compute_stats(adata)

# # schizophrenia
# scz_gene_list=df_gs[df_gs['TRAIT']=='PASS_Schizophrenia_Pardinas2018']['GENESET'].values[0].split(",")
# df_res = scdrs.method.score_cell(adata, scz_gene_list)
# df_res.to_csv("/project2/gilad/jpopp/ebQTL/results/scDRS/scz_results.csv")

# IBD
ibd_gene_list=df_gs[df_gs['TRAIT']=='PASS_IBD_deLange2017']['GENESET'].values[0].split(",")
df_res = scdrs.method.score_cell(adata, ibd_gene_list)
df_res.to_csv("/project2/gilad/jpopp/ebQTL/results/scDRS/ibd_results.csv")

# hypertension
htn_gene_list=df_gs[df_gs['TRAIT']=='UKB_460K.disease_HYPERTENSION_DIAGNOSED']['GENESET'].values[0].split(",")
df_res = scdrs.method.score_cell(adata, htn_gene_list)
df_res.to_csv("/project2/gilad/jpopp/ebQTL/results/scDRS/htn_results.csv")