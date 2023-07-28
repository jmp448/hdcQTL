import pandas as pd
import numpy as np
import tensorqtl
from tensorqtl import genotypeio, cis, trans
from sklearn.preprocessing import scale, quantile_transform

plink_prefix_path = snakemake.params['plink_prefix']
covariates_file = snakemake.input['cov']
expression_bed = snakemake.input['exp']

geno_loc = snakemake.output['genotypes_filtered']
exp_loc = snakemake.output['exp_filtered']
cov_loc = snakemake.output['cov_filtered']

pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
genotype_df_normalized = pd.DataFrame(scale(genotype_df, axis=1), index=genotype_df.index, columns=genotype_df.columns)

covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)

donors_inc = genotype_df_normalized.columns.intersection(covariates_df.columns)
genotype_df_normalized = genotype_df_normalized.loc[:, donors_inc]
covariates_df = covariates_df.loc[:, donors_inc]

phenotype_df, _ = tensorqtl.read_phenotype_bed(expression_bed)
phenotype_df = phenotype_df.loc[:, donors_inc]
phenotype_df_normalized = pd.DataFrame(quantile_transform(phenotype_df, axis=1, output_distribution='normal'),
                                       index=phenotype_df.index, columns=phenotype_df.columns)

genotype_df_normalized.to_csv(geno_loc, sep="\t")
covariates_df.to_csv(cov_loc, sep="\t")
phenotype_df_normalized.to_csv(exp_loc, sep="\t")