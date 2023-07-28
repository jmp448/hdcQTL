import pandas as pd
import torch
import tensorqtl
import numpy as np
from tensorqtl import genotypeio, cis, trans
from sklearn.preprocessing import quantile_transform


print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

# define paths to data
expression_bed = snakemake.input['exp']
covariates_file = snakemake.input['cov']
donor_sample_loc = snakemake.input['donor_sample']
harmonized_tests_loc = snakemake.input['all_tests']
plink_prefix_path = snakemake.params['plink_prefix']
npcs = snakemake.wildcards['npcs']
output_loc = snakemake.params['output_prefix']

# load the tissue-specific donor sample
with open(donor_sample_loc, "r") as file:
    donor_sample = file.readlines()
donor_sample = [s.strip() for s in donor_sample]

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
ix = [pr.sample_ids.index(i) for i in donor_sample]
pr.fam = pr.fam.loc[ix]
pr.bed = pr.bed[:,ix]
pr.sample_ids = pr.fam['iid'].tolist()
genotype_df = pr.load_genotypes()
genotype_df = genotype_df.loc[:, donor_sample]
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# load genes included in the EB analysis
all_tests = pd.read_csv(harmonized_tests_loc, sep="\t")
eb_genes = np.unique(all_tests['phenotype_id_ensg'])

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
phenotype_df.index = phenotype_df.index.str.split('.').str[0]
phenotype_pos_df.index = phenotype_df.index.str.split('.').str[0]

covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0, nrows=int(npcs)+5).T

# filter covariates to sampled donors
covariates_df = covariates_df.loc[donor_sample, :]

# filter phenotypes to sampled donors, overlap genes
phenotype_df = phenotype_df.loc[phenotype_df.index.intersection(eb_genes), donor_sample]
phenotype_pos_df = phenotype_pos_df.loc[phenotype_df.index.intersection(eb_genes), :]

# redo standardization of gene expression
phenotype_df_norm = quantile_transform(phenotype_df, axis=1, output_distribution='normal')
phenotype_df_norm = pd.DataFrame(phenotype_df_norm, index=phenotype_df.index, columns=phenotype_df.columns)

cis.map_nominal(genotype_df, variant_df, 
                phenotype_df_norm, phenotype_pos_df, 
                output_loc, covariates_df,
                maf_threshold=0.1,
                window=50000)
