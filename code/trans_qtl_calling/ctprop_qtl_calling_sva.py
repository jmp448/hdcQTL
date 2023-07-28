import pandas as pd
import torch
import tensorqtl
import numpy as np
from tensorqtl import genotypeio, cis, trans
from sklearn import preprocessing
from sklearn.preprocessing import quantile_transform
import seaborn as sns

genotypes_loc = snakemake.input['genotypes']
proportions_loc = snakemake.input['ctprops']
covariates_loc = snakemake.input['cov']
output_loc = snakemake.output['output_loc']

print(genotypes_loc)
print(proportions_loc)
print(covariates_loc)
print(output_loc)

# Load genotypes
genotype_df = pd.read_csv(genotypes_loc, sep="\t", index_col='snp')

# Load covariates
covariates_df = pd.read_csv(covariates_loc, sep='\t', index_col=0)

# Filter to the donors for whom we have expression + genotype data
donors_inc = genotype_df.columns.intersection(covariates_df.columns)
print(donors_inc)

# Load phenotypes, filter to the donors included and cell types with median score of at least 0.1
phenotype_df = pd.read_csv(proportions_loc, sep="\t").set_index('cell_type')
phenotype_df = phenotype_df.loc[phenotype_df.median(axis=1) > 0.1, donors_inc]

if phenotype_df.shape[0] == 0: # Skip QTL calling in tissues w/o any of these cell types
    pd.DataFrame(columns=['variant_id', 'phenotype_id', 'pval', 'b', 'b_se', 'af']).to_csv(output_loc, index=False)
else:
    # Quantile normalize the phenotypes
    phenotype_df_norm = quantile_transform(phenotype_df, axis=1, output_distribution='normal')
    phenotype_df_norm = pd.DataFrame(phenotype_df_norm, index=phenotype_df.index, columns=phenotype_df.columns)

    # Use variant-specific covariates
    qtl_all = []
    for s in covariates_df.index.unique().tolist():
        geno_df = genotype_df.loc[genotype_df.index == s, :]
        cov_df = covariates_df.loc[covariates_df.index == s, :].set_index('ID').T
        qtl_all.append(trans.map_trans(geno_df, phenotype_df_norm, cov_df,
                                       pval_threshold=1, maf_threshold=0))
    trans_df = pd.concat(qtl_all)

    # Perform cell type proportion QTL calling
    trans_df.to_csv(output_loc, index=False)
