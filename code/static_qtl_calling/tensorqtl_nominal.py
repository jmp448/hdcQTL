import pandas as pd
import torch
import tensorqtl
import numpy as np
import seaborn as sns
from tensorqtl import genotypeio, cis, trans
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

# define paths to data
expression_bed = snakemake.input['exp']
covariates_file = snakemake.input['cov']
plink_prefix_path = snakemake.params['plink_prefix']
output_loc = snakemake.params['output_prefix']

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

cis.map_nominal(genotype_df, variant_df, 
                phenotype_df, phenotype_pos_df, 
                output_loc, covariates_df,
                window=50000)
