import pandas as pd
import torch
import tensorqtl
import numpy as np
from tensorqtl import genotypeio, cis, trans
from sklearn import preprocessing
import seaborn as sns

plink_prefix_path = snakemake.params['plink_prefix']
expression_bed = snakemake.input['exp']
covariates_file = snakemake.input['cov']
output_loc = snakemake.output['output_loc']

# Load genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# Load phenotypes
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)

# Reformat gene names to omit version number
phenotype_df.index = phenotype_df.index.str.split('.').str.get(0)
phenotype_pos_df.index = phenotype_pos_df.index.str.split('.').str.get(0)

# Load covariates
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# Filter to the donors for whom we have expression + genotype data
donors_inc = genotype_df.columns.intersection(covariates_df.index)
genotype_df = genotype_df.loc[:, donors_inc]
covariates_df = covariates_df.loc[donors_inc, :]
phenotype_df = phenotype_df.loc[:, donors_inc]

# assert np.all(phenotype_df.columns==covariates_df.index)
# assert covariates_df.index.isin(genotype_df.columns).all()

# Perform trans eQTL calling
trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df,
                           pval_threshold=1, maf_threshold=0.05,
                           batch_size=20000)

# Filter to tests where variant & chromosome aren't on the same gene
trans_filtered = trans.filter_cis(trans_df, phenotype_pos_df.T.to_dict(), variant_df, window=np.Inf)

trans_filtered.to_csv(output_loc)


