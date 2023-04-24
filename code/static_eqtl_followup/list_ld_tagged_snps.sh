#!/bin/bash

module load vcftools

genotypes="$1"
snps="$2"
prefix="$3"

vcftools --vcf $genotypes --geno-r2-positions $snps \
          --ld-window-bp 50000 \
          --min-r2 0.8 --out $prefix
          
# vcftools --gzvcf "data/genotypes/human.YRI.hg38.all.AF.gencode.vcf.gz" --keep "data/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/individuals.tsv" \
#           --geno-r2-positions "results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/matrixeqtl.cis_qtl_pairs.tophits.pos" \
#           --ld-window-bp 100 --min-r2 0.8s \
#           --out "results/static/ebqtl_ipsc/pseudobulk_tmm/basic/IPSC/8pcs/yri_tagging"