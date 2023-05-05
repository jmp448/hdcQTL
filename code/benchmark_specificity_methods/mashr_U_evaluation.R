library(tidyverse)
library(mashr)

# Load info from snakemake
mash_data_loc <- snakemake@input[[1]]
output_loc <- snakemake@output[[1]]

# Load info locally
load(mash_data_loc)

data.train <- mash_set_data(Bhat.train, Shat.train)
data.test <- mash_set_data(Bhat.test, Shat.test)
data.strong <- mash_set_data(Bhat.strong, Shat.strong)
data.strong.ez <- mash_set_data(Bhat.strong, Shat.strong, alpha=1)

# Collect all potential covariance matrices
U.c.all = cov_canonical(data.train)
U.c.equal_effects <- list("equal_effects"=U.c.all$equal_effects)
U.c.equal_effects_identity <- list("equal_effects"=U.c.all$equal_effects, "identity"=U.c.all$identity)

U.pca.default = cov_pca(data.strong, 5)
U.pca.ez = cov_pca(data.strong.ez, 5)
U.ed.default = cov_ed(data.strong, U.pca.default)
U.flash = cov_flash(data.strong)

U.pca_no_combo = U.pca
U.pca_no_combo[['tPCA']] <- NULL
U.flash_no_combo = U.flash
U.flash_no_combo[['tFLASH_default']] <- NULL

# Baseline round 1 - canonical only
m.equal_effects_only = mash(data.train, U.c.equal_effects, outputlevel = 4)
m.equal_effects_identity =  mash(data.train, U.c.equal_effects_identity, outputlevel = 4)
m.canonical =  mash(data.train, U.c.all, outputlevel = 4)

canonical_models <- list("equal_effects"=m.equal_effects_only,
                         "equal_effects_identity"=m.equal_effects_identity,
                         "canonical"=m.canonical)

canonical_comp <- tibble("equal_effects"=mash_compute_loglik(m.equal_effects_only, data.test),
                         "equal_effects_identity"=mash_compute_loglik(m.equal_effects_identity, data.test),
                         "canonical_all"=mash_compute_loglik(m.canonical, data.test))

# Baseline round 2 - incorporate some data-driven covariance matrices
m.default = mash(data.train, Ulist = c(U.ed.default, U.c.all), outputlevel = 4)
m.pca = mash(data.train, Ulist = c(U.pca.default, U.c.all), outputlevel = 4)
m.pca.ez = mash(data.train, Ulist = c(U.pca.ez, U.c.all), outputlevel = 4)
m.pca_no_combo = mash(data.train, Ulist = c(U.pca_no_combo, U.c.all), outputlevel = 4)
m.flash = mash(data.train, Ulist = c(U.flash, U.c.all), outputlevel = 4)
m.flash_no_combo = mash(data.train, Ulist = c(U.flash_no_combo, U.c.all), outputlevel = 4)

datadriven_models <- list("ed"=m.default,
                          "pca"=m.pca,
                          "pca_ez"=m.pca.ez,
                          "pca_nocombo"=m.pca_no_combo,
                          "flash"=m.flash,
                          "flash_nocombo"=m.flash_no_combo)

datadriven_comp <- tibble("ed"=mash_compute_loglik(m.default, data.test),
                         "pca"=mash_compute_loglik(m.pca, data.test),
                         "pca_ez"=mash_compute_loglik(m.pca.ez, data.test),
                         "pca_nocombo"=mash_compute_loglik(m.pca_no_combo, data.test),
                         "flash"=mash_compute_loglik(m.flash, data.test),
                         "flash_nocombo"=mash_compute_loglik(m.flash_no_combo, data.test))

# Save outputs - train/test data, fitted mash models, 
save(canonical_models, canonical_comp,
     datadriven_models, datadriven_comp,
     file=output_loc)