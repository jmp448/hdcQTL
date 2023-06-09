###
# This file computes MAF using vcftools
# Allows zero missingness, filters to a set of individuals stored in the fourth input file, filters to MAF 0.1
# Also 'thins' to remove adjacent SNPs (this actually just removes one SNP that was duplicated because it has two rsIDs)
# Also filters to autosomal variants (spelled out chromosomes 1-22 because the VCF file didn't just have chrX and chrY in addition, but weird other stuff)
###

module load vcftools

genotypes="$1"
inds="$2"
prefix="$3"
positions="$4"

vcftools --gzvcf $genotypes --out $prefix --keep $inds --positions $positions --freq