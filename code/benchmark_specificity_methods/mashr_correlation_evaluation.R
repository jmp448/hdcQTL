library(tidyverse)
library(mashr)
set.seed(1234)

# Load info from snakemake
mash_data_loc <- snakemake@input[[1]]
output_loc <- snakemake@output[[1]]

# Load info locally
load(mash_data_loc)

# Baseline - assume no correlation between contexts
data.train <- mash_set_data(Bhat.train, Shat.train)
data.test <- mash_set_data(Bhat.test, Shat.test)

# Estimate correlation using the basic approach
V.simple = estimate_null_correlation_simple(data.train)
data.train.Vsimple = mash_update_data(data.train, V=V.simple)
data.test.Vsimple = mash_update_data(data.test, V=V.simple)

U.c = cov_canonical(data.train)

# Estimate mash parameters using basic approach (with and without simple correlation estimate)
m.orig = mash(data.train, U.c)
m.Vsimple = mash(data.train.Vsimple, U.c)

# Estimate correlation & mash parameters using EM approach
V.em.full = mash_estimate_corr_em(data.train, U.c, max_iter=10, details = TRUE)
V.em = V.em.full$V
m.Vem = V.em.full$mash.model
data.test.Vem = mash_update_data(data.test, V=V.em)

# Compare model fits on held-out test data
m.orig.test = mash(data.test, g=get_fitted_g(m.orig), fixg=TRUE)
m.Vsimple.test = mash(data.test.Vsimple, g=get_fitted_g(m.Vsimple), fixg=TRUE)
m.Vem.test = mash(data.test.Vem, g=get_fitted_g(m.Vem), fixg=TRUE)

# Save outputs - train/test data, fitted mash models, 
save(m.orig.test, m.Vsimple.test, m.Vem.test,
     file=output_loc)