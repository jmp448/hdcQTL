#!/bin/bash
#SBATCH --time=5:0:0
#SBATCH --partition=gilad
#SBATCH --job-name=scran
#SBATCH --output=logs/scran.out
#SBATCH --error=logs/scran.err
#SBATCH --mem=500G

module load R/4.0.0
Rscript ../code/normalization_scran.R