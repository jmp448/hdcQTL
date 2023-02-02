###
# Reformats some of the output from vcftools for easier input to Matrix eQTL
# 
# We add ID's to the snp locs file (for now, just using chr1--1 for example, since some variants don't have rsID)
# We also remove the stupid indices from the dosage-encoded genotype file, add the donor name, and add snp ID for the genotype file
###

genotypes="$1"
individuals="$2"
snp_locs="$3"
temploc="$4"
snp_locs_out="$5"
genotypes_out="$6"

### add snp ids (just positions - some variants lack rsID) to snp locs file
sed 's/\t/_/g' $snp_locs > $temploc.snpid.tmp
paste $temploc.snpid.tmp $snp_locs > $snp_locs_out

### add snp ids and individual names to genotype file
# convert individual names to header
bash code/mashr/genotype_transpose.sh $individuals $temploc.genotypes.tmp 

# add genotype values below header
tail -n +2 $genotypes >> $temploc.genotypes.tmp 

# add snp ids to genotypes file
echo snp_id > $temploc.snpid.header.tmp
cat $temploc.snpid.tmp >> $temploc.snpid.header.tmp
paste $temploc.snpid.header.tmp $temploc.genotypes.tmp > $genotypes_out

rm $temploc.*

