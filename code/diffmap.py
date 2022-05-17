import scanpy as sc
import pandas as pd
import numpy as np
import scvi
from glob import glob

adata = sc.read_h5ad("/project2/gilad/katie/ebQTL/highpass_combinedFiles/102andPilot_MetaAdded_QCadded_filtered_noNorm_5000VarFeatNoBatchKey_LinearSCVI022222.h5ad")

sc.tl.diffmap(adata)

adata.write_h5ad("/project2/gilad/jpopp/ebQTL/data/single_cell_objects/highpass/diffmapped.h5ad")
