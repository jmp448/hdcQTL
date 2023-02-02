###
# Use vcftools to convert vcf format to a dosage encoding (0/1/2)
###

module load vcftools

genotypes="$1"
prefix="$2"

vcftools --gzvcf $genotypes --out $prefix --012