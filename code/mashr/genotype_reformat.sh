###
# Reformats some of the output from vcftools for easier input to Matrix eQTL
# 
# We add ID's to the snp locs file (for now, just using chr1--1 for example, since some variants don't have rsID)
# We also remove the stupid indices from the dosage-encoded genotype file, add the donor name, and add snp ID for the genotype file
###

genotypes="$1"
individuals="$2"
snp_locs="$3"
agg=`echo $snp_locs | cut -d/ -f3`
annot=`echo $snp_locs | cut -d/ -f5`
type=`echo $snp_locs | cut -d/ -f6`


snp_locs_out="$4"
genotypes_out="$5"

### add snp ids (just positions - some variants lack rsID) to snp locs file
sed 's/\t/_/g' $snp_locs > snpid.$type.$annot.$agg.tmp
paste snpid.$type.$annot.$agg.tmp $snp_locs > $snp_locs_out

### add snp ids and individual names to genotype file
# convert individual names to header
bash code/mashr/genotype_transpose.sh $individuals genotypes.$type.$annot.$agg.tmp 

# add genotype values below header
tail -n +2 $genotypes >> genotypes.$type.$annot.$agg.tmp 

# add snp ids to genotypes file
echo snp_id > snpid.header.$type.$annot.$agg.tmp
cat snpid.$type.$annot.$agg.tmp >> snpid.header.$type.$annot.$agg.tmp
paste snpid.header.$type.$annot.$agg.tmp genotypes.$type.$annot.$agg.tmp > $genotypes_out

rm *.$type.$annot.$agg.tmp

