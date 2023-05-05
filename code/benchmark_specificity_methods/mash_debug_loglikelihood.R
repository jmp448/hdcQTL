library(mashr)
library(mvtnorm)
library(tidyverse)

mash_data_loc <- "results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/tmmtype/8pcs/tensorqtlgenonorm_mash_inputs.train_test_split.Rdata"
load(mash_data_loc)

# Baseline - assume no correlation between contexts
data.train <- mash_set_data(Bhat.train, Shat.train)
data.test <- mash_set_data(Bhat.test, Shat.test)

U.c = cov_canonical(data.train)

# Estimate mash parameters using basic approach 
m.orig = mash(data.train, U.c, outputlevel=4)

# The log likelihood of individual tests computed here is not valid (should be (-inf, 0] but mostly positive)
hist(m.orig$vloglik)

# Are the relative likelihoods used to compute that sum valid?
llik_matrix <- m.orig$lm$loglik_matrix + m.orig$lm$lfactors
hist(llik_matrix)

# No, these are not valid log likelihoods

# Is this an issue with the scaling procedure, or the original calculation of the log likelihood?
original_m_llik <- mashr:::calc_lik_matrix(data.train, U.c, log=T)
hist(original_m_llik)

# These are also not valid log likelihoods, so it's an issue with the likelihood calculation
# Looking at that likelihood for just a few tests to hone in on the issue
res <- t(sapply(1:100, function(j) mashr:::calc_lik_vector(data.train$Bhat[j,],mashr:::get_cov(data.train,j),
                                            U.c,log=T)))

# Actually, just in the very first test we get an invalid log likelihood vector
test1 <- mashr:::calc_lik_vector(data.train$Bhat[1,],mashr:::get_cov(data.train,1),U.c,log=T)

# Why is this?
## The V matrix looks right
all.equal(diag(mashr:::get_cov(data.train,1)), unname(data.train$Shat[1,]^2))

## The singleton covariance matrices for some reason do not return a valid log likelihood...
dmvnorm(data.train$Bhat[1,], sigma = U.c[[1]] + diag(unname(data.train$Shat[1,]^2)),log = T)
dmvnorm(data.train$Bhat[1,], sigma = U.c[[2]] + diag(unname(data.train$Shat[1,]^2)),log = T)

# This has to do with the way likelihood is defined - dmvnorm offers only a relative likelihood
U.c.std <- U.c
for (i in length(U.c)) {
  U.c.std[[i]] <- U.c[[i]] / max(diag(U.c[[1]]))
}

dmvnorm(data.train$Bhat[1,], sigma = U.c[[29]] + diag(unname(data.train$Shat[1,]^2)),log = T)

# I think this is a non-issue, the relative log likelihood is valid, it just isn't expected to have the 
# same bounds as an absolute likelihood. But are likelihoods comparable as is? Should they be comparable 
# across U's? Are they meaningfully compared?

# For relative to apply, we want the total probability over all possible values of beta to be similar across
# each of our U matrices
as_tibble(original_m_llik) %>%
  pivot_longer(everything(), names_to="component", values_to="RLL") %>%
  filter(component=="simple_het_1") %>%
  ggplot(aes(x=RLL, fill=component)) +
  geom_density() +
  xlim(-50, 50)

# Instead, the components which span multiple cell types have density spiked low, 
# and the components for a single cell type are spread over MUCH larger values

# Would these be more comparable if they were normalized by the area under the curve from
# (uniformly) zero effect to (uniformly) an effect of one?

# Actually, first why don't I figure out if this is actually causing an issue
### The distribution of likelihoods over all tests differed dramatically between methods
data.train.genonorm <- data.train
data.test.genonorm <- data.test
m.orig.genonorm <- m.orig

mash_data_basic_loc <- "results/benchmark_specificity_methods/eb_cellid/pseudobulk_tmm/tmmtype/8pcs/tensorqtlbasic_mash_inputs.train_test_split.Rdata"
load(mash_data_basic_loc)

data.train.basic <- mash_set_data(Bhat.train, Shat.train)
data.test.basic <- mash_set_data(Bhat.test, Shat.test)

m.orig.basic = mash(data.train.basic, U.c, outputlevel=4)

tibble("genotypes"=rep(c("dosage", "normalized"), each=20000), "vloglik"=c(m.orig.basic$vloglik, m.orig.genonorm$vloglik)) %>%
  ggplot(aes(x=vloglik, fill=genotypes)) + geom_density()

# Does this have anything to do with the prior, or do we see it on just p(\hat{b}|V)?
just_se_llik_basic <- mashr:::calc_lik_matrix(data.train.basic, Ulist=list(matrix(0, nrow=29, ncol=29)), log=T)
just_se_llik_genonorm <- mashr:::calc_lik_matrix(data.train.genonorm, Ulist=list(matrix(0, nrow=29, ncol=29)), log=T)

tibble("genotypes"=rep(c("dosage", "normalized"), each=20000), "vloglik"=c(just_se_llik_basic[,1], just_se_llik_genonorm[,1])) %>%
  ggplot(aes(x=vloglik, fill=genotypes)) + geom_density()

# If the standard errors systematically decreased upon genotype normalization, the relative likelihoods would have 
# systematically increased, making this a meaningless comparison
tibble("genotypes"=rep(c("dosage", "normalized"), each=580000), "se"=c(as.numeric(data.train.basic$Shat)^2, as.numeric(data.train.genonorm$Shat)^2)) %>%
  ggplot(aes(x=se, fill=genotypes)) + geom_density()

# Still though, it really seems like a problem that the relative likelihoods are not 
# comparable across distributions when we're trying to maximize the likelihood jointly
# This can be visualized when looking at the likelihood of all tests for a singleton component vs a joint one
as_tibble(original_m_llik) %>%
  pivot_longer(everything(), names_to="component", values_to="RLL") %>%
  filter(component %in% c("equal_effects", "simple_het_1", "Parietal-and-chief-cells")) %>%
  ggplot(aes(x=RLL, fill=component)) +
  geom_density() +
  xlim(-50, 50) +
  theme_classic(base_size=20) +
  theme(legend.position="none") +
  facet_grid(component ~ ., scales="free")
p0_equal_effects <- dmvnorm(rep(0, 29), sigma = U.c[['equal_effects']] + diag(unname(data.train$Shat[1,]^2)), log = T)
p0_simple_het <- dmvnorm(rep(0, 29), sigma = U.c[['simple_het_1']] + diag(unname(data.train$Shat[1,]^2)), log = T)
p0_parietal <- dmvnorm(rep(0, 29), sigma = U.c[['Parietal-and-chief-cells']] + diag(unname(data.train$Shat[1,]^2)), log = T)

# To see if there's evidence of this bias in our real data, need the mash fit from our data
m <- readRDS("results/static_qtl_calling/eb_cellid/pseudobulk_tmm/basic/8pcs/mash_fitted_model.tophits.rds")
weights <- as_tibble(get_estimated_pi(m), rownames="component") 
rll_0_identity <- tibble("component"=weights$component, "rll"=NA_real_)
rll_0_variable_se <- tibble("component"=weights$component, "rll"=NA_real_)
for (i in seq(2, nrow(rll_0))) {
  rll_0_identity[i, "rll"] <- dmvnorm(rep(0, 29), sigma = m$fitted_g$Ulist[[rll_0$component[i]]] + diag(1, nrow=29, ncol=29), log = T)
  rll_0_variable_se[i, "rll"] <- dmvnorm(rep(0, 29), sigma = m$fitted_g$Ulist[[rll_0$component[i]]] + diag(unname(data.train$Shat[1,]^2)), log = T)
}
inner_join(weights, rll_0_identity, by="component") %>%
  drop_na() %>% # this just drops the null component
  ggplot(aes(x=rll, y=value)) +
  geom_point() +
  xlab("RLL at 0") +
  ylab("Assigned weight")

inner_join(weights, rll_0_variable_se, by="component") %>%
  drop_na() %>% # this just drops the null component
  ggplot(aes(x=rll, y=value)) +
  geom_point() +
  xlab("RLL at 0") +
  ylab("Assigned weight")

## This bias actually may be desired behavior, but I need to understand whether my choice of
## data-driven matrices is reducing their chances, as they may not actually be rank-1 which would hurt them