import pandas as pd
import numpy as np
from pandas_plink import read_plink1_bin

test_eqtl_file = snakemake.input['test_eqtl_file']
genotype_file = snakemake.input['genotype_file']
filtered_test_file = snakemake.output['filtered_eqtl_file']

all_tests = pd.read_csv(test_eqtl_file, sep="\t")
G = read_plink1_bin(genotype_file)

maf01snps = set(G.snp.values)
all_tests = all_tests[all_tests['EB_VARIANT_ID'].isin(maf01snps)]

all_tests.to_csv(filtered_test_file, sep="\t", index=False)
