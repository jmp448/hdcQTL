###
# Write genotype data to BED format for tensorQTL
###

module load plink

genotypes="$1"
inds="$2"
prefix="$3"

plink --make-bed \
    --chr 1-22 \
    --maf 0.1 \
    --geno 0 \
    --bp-space 1 \
    --aec \
    --keep-fam ${inds} \
    --keep-allele-order \
    --vcf ${genotypes} \
    --output-chr chrM \
    --out ${prefix}