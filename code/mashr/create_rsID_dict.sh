#!/bin/sh

module load bcftools

bcftools query -f '%CHROM\t%POS\t%ID\n' human.YRI.hg38.all.AF.gencode.vcf.gz

