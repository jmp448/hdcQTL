#!/bin/bash
#SBATCH --time=6:0:0
#SBATCH --partition=gilad
#SBATCH --job-name=scdrs
#SBATCH --output=logs/scdrs.out
#SBATCH --error=logs/scdrs.err
#SBATCH --mem=300G

module load python
source activate scvi-gpflow
python ../code/scDRS.py

