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
output_loc = snakemake.output[0]
npcs = snakemake.wildcards['npcs']

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0, nrows=int(npcs)).T

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

cis_df = cis.map_cis(genotype_df, variant_df, 
                     phenotype_df, phenotype_pos_df,
                     covariates_df=covariates_df, seed=123456)

cis_df.to_csv(output_loc, sep="\t")
