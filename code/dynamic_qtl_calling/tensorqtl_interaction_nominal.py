import pandas as pd
import torch
import tensorqtl
import numpy as np
import os
from tensorqtl import genotypeio, cis, trans
from sklearn import preprocessing

# define paths to data
expression_bed = snakemake.input['exp']
covariates_file = snakemake.input['covariate_df']
pseudotime_file = snakemake.input['pseudotime']
genotype_file = snakemake.input['genotype_df']
variant_file = snakemake.input['variant_df']
output_loc = snakemake.params['output_prefix']
n_clpcs = int(snakemake.wildcards['n_cl_pcs'])

# Create the directory for the outputs if needed 
if not os.path.exists(output_loc.rsplit('/', 1)[0]):
    os.makedirs(output_loc.rsplit('/', 1)[0])

# Load phenotype files
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)

# Load cell line PC covariates
covariates_df = pd.read_csv(covariates_file, sep="\t", index_col="sample")

# Load genotype files
genotypes_df = pd.read_csv(genotype_file, sep="\t", index_col="snp")
variant_df = pd.read_csv(variant_file, sep="\t", index_col="snp")

# Load pseudotime (interaction) data
pseudotime = pd.read_csv(pseudotime_file, sep='\t').set_index('ind_type')

assert (list(pseudotime.index) == phenotype_df.columns).all()
assert (list(covariates_df.index) == phenotype_df.columns).all()
assert (genotypes_df.columns == phenotype_df.columns).all()

for i in range(1,23):
    cis.map_nominal(genotypes_df, variant_df,
                    phenotype_df.loc[phenotype_pos_df['chr']=='chr'+str(i)],
                    phenotype_pos_df.loc[phenotype_pos_df['chr']=='chr'+str(i)], 
                    prefix=output_loc,
                    interaction_df=pseudotime, 
                    covariates_df=covariates_df, 
                    window=50000,
                    maf_threshold_interaction=0,
                    run_eigenmt=False,
                    write_top=False, 
                    write_stats=True)

### VERSION BEFOORE PREP FILE EXISTED
# # define paths to data
# expression_bed = snakemake.input['exp']
# covariates_file = snakemake.input['cov']
# pseudotime_file = snakemake.input['pseudotime']
# plink_prefix_path = snakemake.params['plink_prefix']
# output_loc = snakemake.params['output_prefix']
# n_clpcs = int(snakemake.wildcards['n_cl_pcs'])
# 
# # Create the directory for the outputs if needed 
# if not os.path.exists(output_loc.rsplit('/', 1)[0]):
#     os.makedirs(output_loc.rsplit('/', 1)[0])
# 
# # Load phenotype files
# phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
# 
# # Load and wrangle cell line PC covariates
# cell_line_pc_df = pd.read_csv(covariates_file, sep='\t', usecols=range(n_clpcs + 1))
# 
# cell_line_pc_long = pd.melt(cell_line_pc_df, id_vars=['ind'], var_name='covariate', value_name='value')
# cell_line_pc_long['covariate'] = cell_line_pc_long['covariate'].apply(lambda x: 'CL' + x)
# cell_line_pc_long = cell_line_pc_long.astype({'covariate': 'string', 'ind': 'string'})
# 
# sample_map = pd.DataFrame({'sample': phenotype_df.columns})
# sample_map['ind'] = sample_map['sample'].str.slice(0, 5)
# 
# covariate_map = pd.merge(cell_line_pc_long, sample_map, on='ind', how='right').drop(columns='ind')
# covariates_df = covariate_map.pivot(index='sample', columns='covariate', values='value')
# 
# # Load and wrangle genotype files
# pr = genotypeio.PlinkReader(plink_prefix_path)
# genotypes = pr.load_genotypes()
# variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
# # variant_df['chrom'] = variant_df['chrom'].apply(lambda x: int(x[3:]))
# 
# genotypes = genotypes.rename(columns=lambda x: x[2:])
# genotypes = pd.DataFrame(preprocessing.scale(genotypes, axis=1), index=genotypes.index, columns=genotypes.columns).reset_index()
# 
# # drop duplicated snp names
# duplicates = genotypes.duplicated(subset=['snp'], keep=False)
# print('Dropping ' + str(sum(duplicates)) + ' SNPs with duplicated IDs')
# genotypes = genotypes[~duplicates]
# variant_df = variant_df.loc[genotypes['snp']]
# 
# genotypes = pd.melt(genotypes, id_vars=['snp'], var_name='ind', value_name='genotype')
# genotypes = genotypes.astype({'ind': 'string'})
# 
# genotypes = pd.merge(genotypes, sample_map, on='ind', how='right').drop(columns='ind')
# genotypes_df = genotypes.pivot(index='snp', columns='sample', values='genotype').reindex(variant_df.index)
# 
# # Load pseudotime (interaction) data
# pseudotime = pd.read_csv(pseudotime_file, sep='\t').set_index('ind_type')
# 
# assert (list(pseudotime.index) == phenotype_df.columns).all()
# assert (list(covariates_df.index) == phenotype_df.columns).all()
# assert (genotypes_df.columns == phenotype_df.columns).all()
# 
# for i in range(1,23):
#     cis.map_nominal(genotypes_df, variant_df,
#                     phenotype_df.loc[phenotype_pos_df['chr']=='chr'+str(i)],
#                     phenotype_pos_df.loc[phenotype_pos_df['chr']=='chr'+str(i)], 
#                     prefix=output_loc,
#                     interaction_df=pseudotime, 
#                     covariates_df=covariates_df, 
#                     window=50000, 
#                     run_eigenmt=True,
#                     write_top=False, 
#                     write_stats=True)
# 
# # cis.map_nominal(genotypes_df, variant_df,
# #                 phenotype_df, phenotype_pos_df,
# #                 prefix=output_loc, 
# #                 interaction_df=pseudotime,
# #                 covariates_df=covariates_df,
# #                 window=50000,
# #                 run_eigenmt=True,
# #                 write_top=False,
# #                 write_stats=True)

"""
### RADHIKA ONLY


# Paths to input and output files
data_path = '/project2/gilad/rjangi1/ebQTL/dynamic_eQTL/data/cm/'
expression_bed = data_path+'expression/TMM_norm_pseudobulk.bed.gz'
covariates_file = data_path+'covariates/cell-line_pcs.tsv'
plink_prefix_path = data_path+'plink_genotypes/plink_genotypes'
genotype_matrix_loc = data_path+''

output_dir = data_path+'test_res' 

# Load first 5 cell line PCs and reformat
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T.iloc[:10]

two = covariates_df.copy(deep=True).add_suffix('_2-mid-early')
three = covariates_df.copy(deep=True).add_suffix('_3-mid')
four = covariates_df.copy(deep=True).add_suffix('_4-mid-late')
five = covariates_df.copy(deep=True).add_suffix('_5-late')
covariates_df = covariates_df.add_suffix('_1-early')

covariates_df = pd.concat([covariates_df, two, three, four, five], 
                  axis = 1)


### Workaround if tensorqtl can't read in genotype_df ###
#genotype_df = pd.read_csv(data_path+'genotypes/genotypes_filtered_plink.tsv',sep='\t').rename(columns={'Unnamed: 0': 'snp'}).set_index('snp')

# Reformat genotype matrix to sample level for tensorqtl input
genotype_df.columns = genotype_df.columns.str.replace("NA","")

two = genotype_df.copy(deep=True).add_suffix('_2-mid-early')
three = genotype_df.copy(deep=True).add_suffix('_3-mid')
four = genotype_df.copy(deep=True).add_suffix('_4-mid-late')
five = genotype_df.copy(deep=True).add_suffix('_5-late')
genotype_df = genotype_df.add_suffix('_1-early')

genotype_df = pd.concat([genotype_df, two, three, four, five], 
                  axis = 1)

# Read in median pseudotime as pandas series for interaction term
med_pt = pd.read_csv(data_path+'interaction_terms/median_pseudotime.tsv', sep='\t').set_index('ind_bin')

# Reindex and reformat all inputs to match phenotype_df
med_pt_reindex = med_pt.reindex(phenotype_df.columns)
genotype_df = genotype_df[phenotype_df.columns]
covariates_df = covariates_df[phenotype_df.columns].T

# Run TensorQTL - window size = 50,000, maf_threshold_interaction = 0.1
for i in range(1,23):
    cis.map_nominal(genotype_df, 
                    variant_df,
                    phenotype_df.loc[phenotype_pos_df['chr']=='chr'+str(i)],
                    phenotype_pos_df.loc[phenotype_pos_df['chr']=='chr'+str(i)], 
                    prefix=str(i),
                    interaction_df=med_pt_reindex, 
                    window=50000, 
                    covariates_df=covariates_df, 
                    run_eigenmt=False,
                    output_dir=output_dir, 
                    write_top=False, 
                    write_stats=True)
# Output should save in /project2/gilad/rjangi1/ebQTL/dynamic_eQTL/data/cm/test_res/
"""
