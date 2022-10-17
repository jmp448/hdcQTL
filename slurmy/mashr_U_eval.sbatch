#!/bin/bash
#SBATCH --partition=gilad
#SBATCH --time=10:00:00
#SBATCH --mem=150G
#SBATCH --output=logs/mashr_U_eval/log.out
#SBATCH --error=logs/mashr_U_eval/log.err
#SBATCH --job-name=Ueval

module load R/4.2.0

Rscript code/mashr/mashr_U_evaluation.R
