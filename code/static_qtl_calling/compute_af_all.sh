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

vcftools --gzvcf $genotypes --out $prefix \
                        --max-missing 1 --keep $inds --maf 0.1 --thin 1 \
                        --chr chr1 --chr chr2 --chr chr3 --chr chr4 --chr chr5 --chr chr6 --chr chr7 --chr chr8 \
                        --chr chr9 --chr chr10 --chr chr11 --chr chr12 --chr chr13 --chr chr13 --chr chr14 --chr chr15 \
                        --chr chr16 --chr chr17 --chr chr18 --chr chr19 --chr chr20 --chr chr21 --chr chr22 \
                        --freq