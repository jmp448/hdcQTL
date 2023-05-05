import pandas as pd
import torch
import tensorqtl
import numpy as np
import os
from tensorqtl import genotypeio, cis, trans
from sklearn import preprocessing

# define paths to data
expression_bed = snakemake.input['exp']
covariates_file = snakemake.input['cov']
plink_prefix_path = snakemake.params['plink_prefix']
covariate_output_file = snakemake.output['covariate_df']
genotype_output_file = snakemake.output['genotype_df']
variant_output_file = snakemake.output['variant_df']
output_loc = snakemake.params['output_prefix']
n_clpcs = int(snakemake.wildcards['n_cl_pcs'])

# Create the directory for the outputs if needed 
if not os.path.exists(output_loc.rsplit('/', 1)[0]):
    os.makedirs(output_loc.rsplit('/', 1)[0])

# Load phenotype files
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)

# Load and wrangle cell line PC covariates
cell_line_pc_df = pd.read_csv(covariates_file, sep='\t', usecols=range(n_clpcs + 1))

cell_line_pc_long = pd.melt(cell_line_pc_df, id_vars=['ind'], var_name='covariate', value_name='value')
cell_line_pc_long['covariate'] = cell_line_pc_long['covariate'].apply(lambda x: 'CL' + x)
cell_line_pc_long = cell_line_pc_long.astype({'covariate': 'string', 'ind': 'string'})

sample_map = pd.DataFrame({'sample': phenotype_df.columns})
sample_map['ind'] = sample_map['sample'].str.slice(0, 5)

covariate_map = pd.merge(cell_line_pc_long, sample_map, on='ind', how='right').drop(columns='ind')
covariates_df = covariate_map.pivot(index='sample', columns='covariate', values='value')

# Load and wrangle genotype files
pr = genotypeio.PlinkReader(plink_prefix_path)
genotypes = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

genotypes = genotypes.rename(columns=lambda x: x[2:])
genotypes_normalized = pd.DataFrame(preprocessing.scale(genotypes, axis=1), index=genotypes.index, columns=genotypes.columns).reset_index()

# drop duplicated snp names
duplicates = genotypes_normalized.duplicated(subset=['snp'], keep=False)
print('Dropping ' + str(sum(duplicates)) + ' SNPs with duplicated IDs')
genotypes_normalized = genotypes_normalized[~duplicates]
variant_df = variant_df.loc[genotypes_normalized['snp']]

genotypes_long = pd.melt(genotypes_normalized, id_vars=['snp'], var_name='ind', value_name='genotype')
genotypes_long = genotypes_long.astype({'ind': 'string'})

genotypes_map = pd.merge(genotypes_long, sample_map, on='ind', how='right').drop(columns='ind')
genotypes_df = genotypes_map.pivot(index='snp', columns='sample', values='genotype').reindex(variant_df.index)

# Save genotypes and covariates to a file 
covariates_df.to_csv(covariate_output_file, sep="\t")
variant_df.to_csv(variant_output_file, sep="\t")
genotypes_df.to_csv(genotype_output_file, sep="\t")
