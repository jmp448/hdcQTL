library(tidyverse)
library(mashr)
library(vroom)
set.seed(1234)

mash_input_data_loc <- snakemake@input[['mash_inputs']]
trained_model_loc <- snakemake@input[['trained_model']]

full_fitted_loc <- snakemake@output[['full_fitted_model']]

load(mash_input_data_loc)
m <- readRDS(trained_model_loc)

m_full = mash(data.full, g=get_fitted_g(m), fixg=TRUE, outputlevel=4)
saveRDS(m_full, full_fitted_loc)