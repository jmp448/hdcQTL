#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=06:00:00
#SBATCH --mem=100G
#SBATCH --output=logs/mashr_corr_eval/log.out
#SBATCH --error=logs/mashr_corr_eval/log.err
#SBATCH --job-name=correval

module load R/4.2.0

Rscript code/mashr/mashr_correlation_evaluation.R
