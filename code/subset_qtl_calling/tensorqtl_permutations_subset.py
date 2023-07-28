import pandas as pd
import torch
import tensorqtl
import numpy as np
import seaborn as sns
import random
from tensorqtl import genotypeio, cis, trans
from sklearn.preprocessing import quantile_transform
print(f'PyTorch {torch.__version__}')
print(f'Pandas {pd.__version__}')

# define paths to data
expression_bed = snakemake.input['exp']
covariates_file = snakemake.input['cov']
harmonized_tests_loc = snakemake.input['all_tests']
plink_prefix_path = snakemake.params['plink_prefix']
npcs = snakemake.wildcards['npcs']
donors_loc = snakemake.output['donor_sample']
cis_df_loc = snakemake.output['cis_df']

# set a tissue-customized random seed
tissue_seed = abs(hash(snakemake.wildcards['tissue'])) % (2**32 - 1)
np.random.seed(tissue_seed)

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
phenotype_df.index = phenotype_df.index.str.split('.').str[0]
phenotype_pos_df.index = phenotype_df.index.str.split('.').str[0]

covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0, nrows=int(npcs)+5).T

# sample 53 donors
donor_sample = np.random.choice(phenotype_df.columns.tolist(), size=53, replace=False)

# PLINK reader for genotypes
pr = genotypeio.PlinkReader(plink_prefix_path)
ix = [pr.sample_ids.index(i) for i in donor_sample]
pr.fam = pr.fam.loc[ix]
pr.bed = pr.bed[:,ix]
pr.sample_ids = pr.fam['iid'].tolist()
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

# filter to genes included in the EB analysis
all_tests = pd.read_csv(harmonized_tests_loc, sep="\t")
eb_genes = np.unique(all_tests['phenotype_id_ensg'])

# double check genotype df order
genotype_df = genotype_df.loc[:, donor_sample]

# filter phenotypes to sampled donors, overlap genes
phenotype_df = phenotype_df.loc[phenotype_df.index.intersection(eb_genes), donor_sample]
phenotype_pos_df = phenotype_pos_df.loc[phenotype_df.index.intersection(eb_genes), :]

# filter covariates to sampled donors
covariates_df = covariates_df.loc[donor_sample, :]

# redo standardization of gene expression
phenotype_df_norm = quantile_transform(phenotype_df, axis=1, output_distribution='normal')
phenotype_df_norm = pd.DataFrame(phenotype_df_norm, index=phenotype_df.index, columns=phenotype_df.columns)

# call qtls
cis_df = cis.map_cis(genotype_df, variant_df, 
                     phenotype_df_norm, phenotype_pos_df,
                     covariates_df=covariates_df, 
                     maf_threshold=0.1,
                     window=50000,
                     seed=123456)

# save qtl outputs and the tissue-customized donor sample
cis_df.to_csv(cis_df_loc, sep="\t")

with open(donors_loc, "w") as file:
    for d in list(donor_sample):
        file.write(str(d) + "\n")
