#!/bin/bash
#SBATCH --time=5:0:0
#SBATCH --partition=gilad
#SBATCH --job-name=cellid_assign
#SBATCH --output=../logs/cellid/assign.out
#SBATCH --error=../logs/cellid/assign.err
#SBATCH --mem=500G

module load R/4.2.0
Rscript ../code/annotation/cellid_assign.R
