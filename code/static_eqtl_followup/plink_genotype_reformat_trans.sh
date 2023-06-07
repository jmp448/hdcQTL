###
# Write genotype data to BED format for tensorQTL
###

module load plink

genotypes="$1"
inds="$2"
trans_candidates="$3"
prefix="$4"

plink --make-bed \
    --maf 0.1 \
    --aec \
    --keep-fam ${inds} \
    --keep-allele-order \
    --extract ${trans_candidates} \
    --vcf ${genotypes} \
    --output-chr chrM \
    --out ${prefix}