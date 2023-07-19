# Before VCF files can be used they need to be compressed using bgzip and indexed with a tabix
# Install tabix: 
    # Conda activate YOUR_ENV 
    # conda install -c bioconda tabix
# Use tabix: tabix -p vcf human.YRI.hg38.all.AF.gencode.vcf.gz

# Reference: https://github.com/single-cell-genetics/limix_qtl/wiki/Inputs#genotype-file

module load plink

all_genotypes="$1"
inds="$2"
prefix="$3"
pruned_bed="$4"

# Make a BED file from a pruned set of variants
plink --make-bed \
    --keep-fam $inds \
    --chr 1-22 \
    --maf 0.1 \
    --hwe 1e-6\
    --geno 0 \
    --bp-space 1 \
    --aec \
    --keep-allele-order \
    --indep-pairwise 5000 50 0.15 \
    --vcf $all_genotypes \
    --output-chr chrM \
    --out $pruned_bed

# Make IBD matrix: 
plink --make-rel square --bfile $pruned_bed --extract ${pruned_bed}.prune.in --out $prefix