###
# Use vcftools to convert vcf format to a dosage encoding (0/1/2)
###

module load vcftools

genotypes="$1"
aggregation="$2"
annotation="$3"
type="$4"

vcftools --gzvcf $genotypes --out data/static/$aggregation/type/$annotation/$type/genotypes_filtered --012