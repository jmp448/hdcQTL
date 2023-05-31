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

cis.map_nominal(genotypes_df, variant_df,
                phenotype_df, phenotype_pos_df,
                prefix=output_loc,
                interaction_df=pseudotime,
                covariates_df=covariates_df,
                window=50000,
                maf_threshold_interaction=0,
                run_eigenmt=True,
                write_top=True,
                write_stats=True)
