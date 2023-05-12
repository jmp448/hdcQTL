#! /usr/bin/env Rscript
#
# Adapted from Peter Carbonetto: prefit_poisson_nmf_purified_pbmc.R
# which Pre-fit a Poisson non-negative factorization to the purified PBMC
# single-cell RNA-seq data of Zheng et al (2017).
#
# Here, pre-fits poisson nmf to data from human embryoid bodies
#
# This script is intended to be run from the command-line shell, with
# options that are processed with the optparse package. For example,
# to fit a rank-4 Poisson non-negative matrix factorization by running
# SCD updatees with extrapolation, 2 threads, and with results saved
# to test.rds, run this command:
#
#   ./prefit_poisson_nmf_purified_pbmc.R -k 4 --nc 2 -o test.rds
#
# Running the script without specifying any options will pre-fit a
# rank-3 Poisson non-negative matrix factorization by running EM
# updates (without extrapolation), without multithreading, and with
# results saved to out.rds.
#

# Load a few packages.
library(Matrix)
library(fastTopics)
library(parallel)

# Get inputs/ outputs from snakemake
ncores <- min(as.integer(snakemake@threads), detectCores())
sce_loc <- snakemake@input[['pseudocells_sce']]
k <- as.integer(snakemake@wildcards[['k']])
fasttopics_fit_loc <- snakemake@output[['fasttopics_fit']]

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
# Load the previously prepared data.
cat("Loading data.\n")
pseudocell_exp <- readRDS(sce_loc)
counts <- counts(pseudocell_exp)
cat(sprintf("Loaded %d x %d counts matrix.\n", nrow(counts), ncol(counts)))

# PRE-FIT POISSON NON-NEGATIVE MATRIX FACTORIZATION
# ---------------------------------------------
cat("Running 400 EM updates without extrapolation to find a good initialization.\n")
timing <- system.time(
  fit0 <- fit_poisson_nmf(counts, k = k, numiter = 400, method = "em",
                         control = list(numiter = 4, nc=ncores)))
cat(sprintf("Computation took %0.2f seconds.\n",timing["elapsed"]))

cat("Running 200 scd updates with extrapolation.\n")
fit <- fit_poisson_nmf(counts, fit0=fit0, numiter = 200, method = "scd",
                       control = list(extrapolate = TRUE, numiter = 4, nc=ncores))


# SAVE RESULTS
# ------------
cat("Saving results.\n")
saveRDS(fit, fasttopics_fit_loc)

# SESSION INFO
# ------------
print(sessionInfo())
