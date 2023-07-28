###
# Write genotype data to BED format for tensorQTL
###

module load plink

genotypes="$1"
prefix="$2"

plink --make-bed \
    --aec \
    --keep-allele-order \
    --vcf ${genotypes} \
    --output-chr chrM \
    --out ${prefix}