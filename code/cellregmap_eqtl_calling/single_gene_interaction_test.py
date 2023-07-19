import os
import sys
import pandas as pd
import xarray as xr
from numpy import ones
from numpy.linalg import cholesky
from pandas_plink import read_plink1_bin
from limix.qc import quantile_gaussianize
from sklearn.preprocessing import StandardScaler 
from cellregmap import run_interaction

# INPUT AND OUTPUT
# Input: gene, test_eqtl_file, sample_mapping_file, genotype_file, kinship_file, phenotype_file, cell_context_file

#covariates_file = snakemake.input['covariate_df']
test_eqtl_file = snakemake.input['test_eqtl_file']
sample_mapping_file = snakemake.input['sample_mapping_file']
genotype_file = snakemake.input['genotype_file']
kinship_file = snakemake.input['kinship_file']
phenotype_file = snakemake.input['exp']
cell_context_file = snakemake.input['cell_contexts']

g = str(snakemake.wildcards['g'])

output_loc = snakemake.output['out']

# ############################################
# ############### Gene file ##################
# ############################################


######################################
#####Filter on specific gene-SNP pairs
all_tests = pd.read_csv(test_eqtl_file, sep="\t")

############################################
########## Sample mapping file #############
############################################

## this file will map pseudocells to donors 
## it will also only include donors we have single-cell data for
sample_mapping = pd.read_csv(sample_mapping_file, sep="\t").sort_values(by=['donor_id', 'pseudocell'])

## donor_id are donor IDs, as found in the genotype matrix (G) and GRM covariance (K)
## cell are cell IDs, as found in the scRNA-seq phenotype vector (y) and cell context covariance (C)

## extract unique individuals
donors = sample_mapping["donor_id"].unique()
print("Number of unique donors: {}".format(len(donors)))

############################################
############# Kinship matrix ###############
############################################

## read in GRM (genotype relationship matrix; kinship matrix)
#kinship_file=input_files_dir+"kinship_file.csv"
K = pd.read_csv(kinship_file, sep = "\t", index_col = 0)
assert all(K.columns == K.index) #symmetric matrix, donors x donors

# Make the dataframe as DataArray
K = xr.DataArray(K.values, dims=["donor_0", "donor_1"], coords={"donor_0": K.columns, "donor_1": K.index})
K = K.sortby("donor_0").sortby("donor_1")

# Find the shared donors between kinship and expression data
donors = sorted(set(list(K.donor_0.values)).intersection(donors))
print("Number of donors after kinship intersection: {}".format(len(donors)))

# subset to relevant donors
K = K.sel(donor_0=donors, donor_1=donors)
assert all(K.donor_0 == donors)
assert all(K.donor_1 == donors)

# Decompose such as K = hK @ hK.T
hK = cholesky(K.values)
hK = xr.DataArray(hK, dims=["sample", "col"], coords={"sample": K.donor_0.values})
assert all(hK.sample.values == K.donor_0.values)

del K

############################################
##### expand from donors to cells ##########

# Expand hK from donors to cells
hK_expanded = hK.sel(sample=sample_mapping["donor_id"].values)
assert all(hK_expanded.sample.values == sample_mapping["donor_id"].values)


#####################################
############ Phenotypes #############
#####################################

# Phenotype (pseudocell expression)
phenotype = xr.open_dataarray(phenotype_file, autoclose=True)
assert all(phenotype.pseudo_cell.values == sample_mapping["pseudocell"].values)
######################################
########## Cell contexts #############
######################################

# cellular environments
# cells by FastTopic loadings (10)
#C_file = input_files_dir+"k10_topic_byCell.csv"
C = pd.read_csv(cell_context_file, index_col = 0, sep = "\t")

# Select only k1 (Axoneme assembly), k2(embryonic skeletal dev), k3 (innate system dev), k4 (autonomic nervous system dev), k6 (neruon fate specification), k9 (mitosis) topics for context
C = C[['k1','k2', 'k3','k4','k6', 'k9']]

C = xr.DataArray(C.values, dims=["pseudocell", "topic"], coords={"pseudocell": C.index.values, "topic": C.columns.values})
assert all(C.pseudocell.values == sample_mapping["pseudocell"].values)

# normalize cell contexts
scaler = StandardScaler()
scaler.fit(C.values)
C.values = scaler.transform(C.values)

######################################
############ Covariates ##############
######################################

cov = sample_mapping[['pseudocell', 'sex']].set_index('pseudocell')
cov['intercept'] = 1

#####################################
############ Genotypes ##############
#####################################

## read in genotype file (plink format)
#plink_file = snakemake_input[1]
G = read_plink1_bin(genotype_file)
# Select snps appearing for that gene
test_snp = all_tests[all_tests['EB_HGNC']==g]['EB_VARIANT_ID'].values[0]

G_expanded = G[:,G['snp']==test_snp].sel(sample=sample_mapping["donor_id"].values)
assert all(hK_expanded.sample.values == G_expanded.sample.values)

print("G_tested shape is {}".format(G_expanded.shape))

######################################
############# Run CRM ################
######################################
#For each gene on specific chromosome

""" Input list: 
    y: n x 1 (only one gene tested at a time)
    W: n x c, where c is the number of fixed effect covariates (e.g., age, sex..)
    C: n x k, where k is the number of contexts to test for interactions
    G: n x s, where s is the number of SNPs to be tested for a given gene
    hK: n x p, where p is the number of individuals, decomposition of the n x n kinship matrix K
"""

# expression for each gene
y = phenotype.sel(gene=g)
y = quantile_gaussianize(y)
y = y.values.reshape(y.shape[0],1)

W = cov.values
C_val = C.values
G_val = G_expanded.values
hK_val = hK_expanded.values

print("Running for gene {}".format(g))

# run interaction test using CellRegMap
pvals = run_interaction(y=y, W=W, G=G_val, E=C_val, hK=hK_val)[0]

pv = pd.DataFrame({"EB_HGNC":g,
                  "EB_VARIANT_ID":test_snp,
                  "pv":pvals})
pv.head()

pv.to_csv(output_loc, sep="\t", header=False, index=False)


# Adapted from https://github.com/annacuomo/CellRegMap_analyses/blob/df01d14de6813d2d7e9f701313af608336b841e5/endodiff/usage/scripts/association_test_for_one_gene.py#L103
