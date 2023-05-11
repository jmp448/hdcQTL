#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import scanpy as sc
import pandas as pd
import anndata as ad
import numpy as np
from scipy.sparse import csr_matrix


# In[ ]:


adata= sc.read_10x_mtx("/project2/gilad/katie/ebQTL/highpass_combinedFiles/HighpassAgg_All102LibsAndPilot_noNorm/outs/count/filtered_feature_bc_matrix/")


# In[ ]:


adata.X = csr_matrix(adata.X)
adata.write(filename= "/project2/gilad/katie/ebQTL/highpass_combinedFiles/102andPilot_NoMeta_noNorm.h5ad")

