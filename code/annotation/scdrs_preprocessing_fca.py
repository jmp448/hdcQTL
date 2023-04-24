import numpy as np
import pandas as pd
import random
import scipy
import scanpy as sc
import scdrs

fca_loc = snakemake.input['fca']
preprocessed_loc = snakemake.output['preprocessed']

fca = sc.read_h5ad(fca_loc)

scdrs.preprocess(fca, adj_prop='celltype')
sc.pp.neighbors(fca, use_rep='MCA')

fca.write_h5ad(preprocessed_loc)