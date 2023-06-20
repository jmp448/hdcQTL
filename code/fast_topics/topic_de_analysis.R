#! /usr/bin/env Rscript
#
# Here, DE analysis to allow for partial membership to topics in the topic model
#
# This script is intended to be run from the command-line shell, with
# options that are processed with the optparse package. 
#
#

# Load a few packages.
library(Matrix)
require(fastTopics)
library(SingleCellExperiment)
library(parallel)


# Get inputs/ outputs from snakemake
ncores <- min(as.integer(snakemake@threads), detectCores())
sce_loc <- snakemake@input[['pseudocells_sce']]
fasttopics_fit_loc <- snakemake@input[['fasttopics_fit']]
outfile <- snakemake@output[['de_analysis']]

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
# Load the previously prepared data.
cat("Loading data.\n")

# Load the count data
pseudocell_exp <- readRDS(sce_loc)
counts <- t(as(assay(pseudocell_exp), "sparseMatrix"))

# Load the fitting data
fit <- readRDS(fasttopics_fit_loc)

cat("Recover the topic model.\n")
fit_multinom <- poisson2multinom(fit)


# DIFFENTIAL EXPRESSION ANALYSIS 
# ---------------------------------------------
cat("DE analysis to allow for partial membership to groups.\n")
timing <- system.time(
  dfa_out <- de_analysis(fit=fit_multinom, X=counts, control=list(ns=20000, nc=ncores))
)

# SAVE RESULTS
# ------------
cat("Saving results.\n")
save(list(fit = fit, fit_multinom=fit_multinom, counts = counts, dfa_out = dfa_out),
     file = outfile)

# SESSION INFO
# ------------
print(sessionInfo())

