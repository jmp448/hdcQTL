import pandas as pd
import torch
import tensorqtl
import numpy as np
from tensorqtl import genotypeio, cis, trans
from sklearn import preprocessing
from sklearn.preprocessing import quantile_transform
import seaborn as sns

plink_prefix_path = snakemake.params['plink_prefix']
proportions_loc = snakemake.input['ctprops']
covariates_file = snakemake.input['cov']
output_loc = snakemake.output['output_loc']

# Load genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# Load covariates
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# Filter to the donors for whom we have expression + genotype data
donors_inc = genotype_df.columns.intersection(covariates_df.index)
genotype_df = genotype_df.loc[:, donors_inc]
covariates_df = covariates_df.loc[donors_inc, :]

# Load phenotypes, filter to the donors included and cell types with median score of at least 0.1
phenotype_df = pd.read_csv(proportions_loc, sep="\t").set_index('cell_type')
phenotype_df = phenotype_df.loc[phenotype_df.median(axis=1) > 0.1, donors_inc]

if phenotype_df.shape[0] == 0: # Skip QTL calling in tissues w/o any of these cell types
    pd.DataFrame(columns=['variant_id', 'phenotype_id', 'pval', 'b', 'b_se', 'af']).to_csv(output_loc, index=False)
else:
    # Quantile normalize the phenotypes
    phenotype_df_norm = quantile_transform(phenotype_df, axis=1, output_distribution='normal')
    phenotype_df_norm = pd.DataFrame(phenotype_df_norm, index=phenotype_df.index, columns=phenotype_df.columns)

    # assert np.all(phenotype_df.columns==covariates_df.index)
    # assert covariates_df.index.isin(genotype_df.columns).all()

    # Perform cell type proportion QTL calling
    trans_df = trans.map_trans(genotype_df, phenotype_df_norm, covariates_df,
                               pval_threshold=1, maf_threshold=0.05,
                               batch_size=20000)
    trans_df.to_csv(output_loc, index=False)
