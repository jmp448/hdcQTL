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
U.c = cov_canonical(data.train)

logliks <- tibble("alpha"=c(0, 0.25, 0.5, 0.75, 1),
                  "loglik"=rep(0, 5))

training_model_list <- list()
testing_model_list <- list()

# Try multiple possible values for alpha
for (i in seq(1, nrow(logliks))) {
  alpha_i = logliks$alpha[i]
  data.train = mash_set_data(Bhat.train, Shat.train, alpha=alpha_i)
  m.train = mash(data.train, U.c)
  
  data.test = mash_set_data(Bhat.test, Shat.test, alpha=alpha_i)
  m.test = mash(data.test, g=get_fitted_g(m.train), fixg=TRUE)
  logliks$loglik[i] = get_loglik(m.test)
  
  training_model_name <- paste0("train_alpha_", alpha_i)
  training_model_list[[training_model_name]] = m.train
  
  test_model_name <- paste0("test_alpha_", alpha_i)
  testing_model_list[[test_model_name]] = m.test
}

# Save outputs
save(training_model_list, testing_model_list, logliks, 
     file=output_loc)