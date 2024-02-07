library(tidyverse)
library(mashr)
library(vroom)
set.seed(1234)

mash_input_data_loc <- snakemake@input[['mash_inputs']]

trained_model_loc <- snakemake@output[['trained_model']]
tophit_fitted_loc <- snakemake@output[['tophits_fitted_model']]

load(mash_input_data_loc)

m = mash(data.random, Ulist = c(U.flash,U.c), outputlevel = 4)
saveRDS(m, trained_model_loc)

m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE, outputlevel=4)
saveRDS(m2, tophit_fitted_loc)