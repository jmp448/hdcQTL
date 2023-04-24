###
# This file lists all individuals from the master VCF file
###

module load vcftools

vcf-query -l $1 > $2