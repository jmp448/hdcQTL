#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=06:00:00
#SBATCH --mem=100G
#SBATCH --output=logs/mashr_alpha_eval/log.out
#SBATCH --error=logs/mashr_alpha_eval/log.err
#SBATCH --job-name=alphaeval

module load R/4.2.0

Rscript code/mashr/mashr_alpha_evaluation.R
