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
library(optparse)
require(fastTopics)

# Process the command-line arguments.
parser <- OptionParser()
parser <- add_option(parser,c("--countfile","-c"),type = "character",default = "")
parser <- add_option(parser,c("--fit","-f"),type = "character",default = "")
parser <- add_option(parser,c("--out","-o"),type="character",default="out.rds")
parser <- add_option(parser,"--nc",type = "integer",default = 1)
out    <- parse_args(parser)
countfile <- out$countfile 
fitfile   <- out$fit
outfile   <- out$out
nc        <- out$nc

rm(parser,out)

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# LOAD DATA
# ---------
# Load the previously prepared data.
cat("Loading data.\n")

# Load the count data
load(countfile)
counts<-counts

# Load the fitting data
fit <- readRDS(fitfile)$fit

cat("Recover the topic model.\n")
fit_multinom <- poisson2multinom(fit)


# DIFFENTIAL EXPRESSION ANALYSIS 
# ---------------------------------------------
cat("DE analysis to allow for partial membership to groups.\n")
timing <- system.time(
  dfa_out <- de_analysis(fit_multinom,counts, control=list(ns=20000, nc=nc)))

# SAVE RESULTS
# ------------
cat("Saving results.\n")
saveRDS(list(fit = fit, fit_multinom=fit_multinom, counts = counts, dfa_out = dfa_out),file = outfile)

# SESSION INFO
# ------------
print(sessionInfo())

